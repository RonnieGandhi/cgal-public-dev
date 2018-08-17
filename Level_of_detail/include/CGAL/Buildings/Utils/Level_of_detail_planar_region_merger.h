#ifndef CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_MERGER_H
#define CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_MERGER_H

// STL includes.
#include <map>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <cassert>
#include <stdlib.h>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>

// New CGAL includes.
#include <CGAL/Buildings/Utils/Level_of_detail_local_mesh_builder.h>
#include <CGAL/Region_growing/Level_of_detail_planar_region_growing.h>
#include <CGAL/Region_growing/Level_of_detail_facets_based_region_growing.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding>
		class Level_of_detail_planar_region_merger {

		public:
			using Kernel   = InputKernel;
            using Building = InputBuilding;

			using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Vector_3   = typename Kernel::Vector_3;
            using Triangle_2 = typename Kernel::Triangle_2;

			using CDT          = typename Building::CDT;
            using Clean_facet  = typename Building::Clean_facet;
			using Clean_facets = typename Building::Clean_facets;

            using Mesh = CGAL::Polyhedron_3<Kernel>;

            using HDS                 = typename Mesh::HalfedgeDS;
            using Mesh_facet_handle   = typename Mesh::Facet_const_handle;
            using Halfedge_handle     = typename Mesh::Halfedge_const_handle;
            using Facet_vertex_handle = typename Mesh::Facet::Vertex_const_handle;

            using Input_region  = std::vector<Mesh_facet_handle>;
            using Input_regions = std::vector<Input_region>;
            using Mesh_facets   = Input_regions;

            using Local_builder         = CGAL::LOD::Local_mesh_builder<Kernel, HDS, Building, Mesh_facet_handle>;
            using Planar_region_growing = CGAL::LOD::Level_of_detail_planar_region_growing<Kernel, Mesh, Mesh_facets>;

			using All_input_regions = std::vector<Input_regions>;

            using Region_facet   = std::vector<Point_3>;
            using Region_facets  = std::vector<Region_facet>;
            using Output_regions = std::vector<Region_facets>;

            typename Kernel::Compute_squared_length_3         squared_length_3;
            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Color = CGAL::Color;

            using Edge           = typename CDT::Edge;
            using Face_handle    = typename CDT::Face_handle;
            using Vertex_handle  = typename CDT::Vertex_handle;
            using Faces_iterator = typename CDT::Finite_faces_iterator;

            using Vertex_handles = std::vector< std::vector<Vertex_handle> >;
            
            using Final_constraint  = std::pair<Vertex_handle, Vertex_handle>;
            using Final_constraints = std::vector<Final_constraint>;

            using Facets_based_region_growing = CGAL::LOD::Level_of_detail_facets_based_region_growing<Kernel>;

			Level_of_detail_planar_region_merger() :
            m_tolerance(FT(1) / FT(100000)),
            m_use_all_facets(false),
            m_use_old_rg(false),
            m_max_num_iters(100),
            m_area_tolerance(FT(1) / FT(100)),
            m_default_height(-FT(100000000000000)) { 

                srand(time(NULL));
            }

			void merge(Clean_facets &clean_facets) const {

                Output_regions output_regions;
                
                if (m_use_old_rg) create_data_old(clean_facets, output_regions);
                else create_data_new(clean_facets, output_regions);

                if (output_regions.size() == 0) return;
                merge_facets(output_regions, clean_facets);
            }

            void use_all_facets(const bool new_state) {
                m_use_all_facets = new_state;
            }

		private:
            const FT m_tolerance;
            
                  bool m_use_all_facets;
            const bool m_use_old_rg;

            const size_t m_max_num_iters;
            const FT m_area_tolerance;

            const FT m_default_height;
			
            void create_data_old(const Clean_facets &clean_facets, Output_regions &output_regions) const {

                Mesh mesh; Input_regions input_regions;
                create_mesh_and_input_regions(clean_facets, mesh, input_regions);
                create_output_regions(input_regions, output_regions);
            }

            void create_data_new(const Clean_facets &clean_facets, Output_regions &output_regions) const {

                Facets_based_region_growing region_growing(clean_facets);
                region_growing.create_regions(output_regions);
            }

            void create_mesh_and_input_regions(const Clean_facets &clean_facets, Mesh &mesh, Input_regions &input_regions) const {

                mesh.clear();

                Local_builder builder(clean_facets);
                mesh.delegate(builder);

                Mesh_facets mesh_facets;
                builder.get_facets(mesh_facets);

                Planar_region_growing planar_region_growing(mesh_facets);    
                planar_region_growing.use_global_conditions(false);

				All_input_regions all_input_regions;
				planar_region_growing.find_regions(all_input_regions);

                CGAL_postcondition(all_input_regions.size() == 1);
                
                input_regions.clear();
                input_regions = all_input_regions[0];
            }

            void create_output_regions(const Input_regions &input_regions, Output_regions &output_regions) const {
                
                output_regions.clear();
                output_regions.resize(input_regions.size());

                for (size_t i = 0; i < input_regions.size(); ++i)
                    add_output_region(input_regions[i], output_regions[i]);
            }

            void add_output_region(const Input_region &input_region, Region_facets &region_facets) const {
                
                CGAL_precondition(input_region.size() > 0);
                region_facets.clear();

                for (size_t i = 0; i < input_region.size(); ++i)
                    add_output_region_facet(input_region[i], region_facets);
            }

            void add_output_region_facet(const Mesh_facet_handle &fh, Region_facets &region_facets) const {

                Region_facet region_facet;

                Halfedge_handle curr = fh->halfedge();
                Halfedge_handle end  = curr;

                Facet_vertex_handle vh = curr->vertex();
                region_facet.push_back(vh->point());

                curr = curr->next();
                while (curr != end) {

                    vh = curr->vertex();
                    region_facet.push_back(vh->point());

                    curr = curr->next();
                };

                region_facets.push_back(region_facet);
            }

            void merge_facets(Output_regions &output_regions, Clean_facets &clean_facets) const {

                const size_t out_size = output_regions.size();

                clean_facets.clear();
                for (size_t i = 0; i < output_regions.size(); ++i) {

                    // if (i == 0)
                    merge_region_facets(output_regions[i], clean_facets);

                    /*if (!m_use_all_facets) {
                        if (output_regions[i].size() == 2)
                            merge_region_facets(output_regions[i], clean_facets);
                    }
                    else merge_region_facets(output_regions[i], clean_facets); */
                }

                const size_t clean_size = clean_facets.size();
                if (out_size == 0 || out_size > clean_size) {

                    std::cout << "out: "   << out_size   << std::endl;
                    std::cout << "clean: " << clean_size << std::endl << std::endl;
                }
            }

            bool merge_region_facets(Region_facets &region_facets, Clean_facets &clean_facets) const {
                
                const Color color = generate_random_color();

                Vector_3 source_normal;
                const bool success = compute_source_normal(region_facets, source_normal);

                if (!success) {
                 
                    for (size_t i = 0; i < region_facets.size(); ++i)
                        clean_facets.push_back(std::make_pair(region_facets[i], color));
                    return;
                }

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                if (source_normal == -target_normal) source_normal = target_normal;

                FT angle; Vector_3 axis;
                compute_angle_and_axis(source_normal, target_normal, angle, axis);

                Point_3 b;
                compute_barycentre(region_facets, b);

                rotate_region_facets(region_facets, angle, axis, b);

                CDT cdt;
                triangulate_region_facets(region_facets, cdt);
                if (cdt.number_of_faces() != 0) get_back_region_facets(cdt, region_facets);

                rotate_region_facets(region_facets, -angle, axis, b, true);
                
                for (size_t i = 0; i < region_facets.size(); ++i)
                    clean_facets.push_back(std::make_pair(region_facets[i], color));
            }

            void rotate_region_facets(Region_facets &region_facets, const FT angle, const Vector_3 &axis, const Point_3 &b, const bool show = false) const {

                // Region_facets tmp = region_facets;

                for (size_t i = 0; i < region_facets.size(); ++i)
                    rotate_region_facet(region_facets[i], angle, axis, b);

                /*
                if (show) {
                    for (size_t i = 0; i < region_facets.size(); ++i) {

                        for (size_t j = 0; j < region_facets[i].size(); ++j) {
                            if (region_facets[i][j].z() < 250.0) {
                             
                                std::cout << "old: " << tmp[i][j] << std::endl;
                                std::cout << "new: " << region_facets[i][j] << std::endl;
                                std::cout << "here: " << angle << "; " << axis << "; " << b << std::endl;
                            }
                        }
                    }
                } */
            }

            void rotate_region_facet(Region_facet &region_facet, const FT angle, const Vector_3 &axis, const Point_3 &b) const {

                if (angle != FT(0)) 
                    rotate_facet(b, angle, axis, region_facet);
            }

            bool compute_source_normal(const Region_facets &region_facets, Vector_3 &normal) const {
                
                CGAL_precondition(region_facets.size() > 0);
                if (region_facets.size() == 0) return false;

                Vector_3 tmp_normal; 
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < region_facets.size(); ++i) {
                    const bool success = compute_facet_normal(region_facets[i], tmp_normal);

                    if (!success) tmp_normal = Vector_3(FT(0), FT(0), FT(0));

                    x += tmp_normal.x();
                    y += tmp_normal.y();
                    z += tmp_normal.z();
                }

                x /= static_cast<FT>(region_facets.size());
                y /= static_cast<FT>(region_facets.size());
                z /= static_cast<FT>(region_facets.size());

                normal = Vector_3(x, y, z);
                const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));

                if (are_equal_points(normal, zero)) return false;
                // std::cout << normal << " ::: ";

                normalize(normal);
                // std::cout << normal << std::endl;
                return true;
            }

            bool compute_facet_normal(const Region_facet &region_facet, Vector_3 &normal) const {
                
                CGAL_precondition(region_facet.size() >= 3);
                if (region_facet.size() < 3) {
                 
                    return false;

                    std::cout << "error: facet has no vertices" << std::endl;
                    exit(0);
                }

                const bool success = compute_cross_product(region_facet, normal);
                if (success) {
                    
                    normalize(normal);
                    return true;
                }
                return false;
            }

            bool compute_cross_product(const Region_facet &region_facet, Vector_3 &normal) const {

                const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
                for (size_t i = 0; i < region_facet.size(); ++i) {

                    const size_t ip  = (i + 1) % region_facet.size();
                    const size_t ipp = (i + 2) % region_facet.size();

                    const Point_3 &p1 = region_facet[i];
                    const Point_3 &p2 = region_facet[ip];
                    const Point_3 &p3 = region_facet[ipp];

                    const Vector_3 v1 = Vector_3(p2, p1);
                    const Vector_3 v2 = Vector_3(p2, p3);

                    normal = cross_product_3(v1, v2);
                    if (!are_equal_points(normal, zero)) {
                     
                        // std::cout << "norm: " << normal << "; ";
                        return true;
                    }
                }
                return false;

                std::cout << "error: normal not found" << std::endl;
                exit(1);
            }

            void compute_target_normal(Vector_3 &normal) const {
                normal = Vector_3(FT(0), FT(0), FT(1));
            }

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {

				const auto cross = cross_product_3(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const FT dot     = dot_product_3(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                if (angle == FT(0)) return;

                CGAL_precondition(length != FT(0));
                if (length == FT(0)) {

                    std::cout << "error: length = 0" << std::endl;
                    exit(0);
                }
                axis = cross / length;

                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle > half_pi) {
                    
                    angle = static_cast<FT>(CGAL_PI) - angle;
                    axis = -axis;
                }
			}

            void compute_barycentre(const Region_facets &region_facets, Point_3 &b) const {

                CGAL_precondition(region_facets.size() > 0);

                Point_3 tmp_b; FT x = FT(0), y = FT(0), z = FT(0);
                for (size_t i = 0; i < region_facets.size(); ++i) {
                    
                    const Region_facet &region_facet = region_facets[i];                    
                    compute_facet_barycentre(region_facet, tmp_b);

                    x += tmp_b.x();
                    y += tmp_b.y();
                    z += tmp_b.z();
                }

                x /= static_cast<FT>(region_facets.size());
                y /= static_cast<FT>(region_facets.size());
                z /= static_cast<FT>(region_facets.size());

                b = Point_3(x, y, z);
            }

            void compute_facet_barycentre(const Region_facet &region_facet, Point_3 &b) const {
                
                CGAL_precondition(region_facet.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < region_facet.size(); ++i) {
                    const Point_3 &p = region_facet[i];

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(region_facet.size());
                y /= static_cast<FT>(region_facet.size());
                z /= static_cast<FT>(region_facet.size());

                b = Point_3(x, y, z);
            }

            void rotate_facet(const Point_3 &b, const FT angle, const Vector_3 &axis, Region_facet &region_facet) const {

                Point_3 q;
                for (size_t i = 0; i < region_facet.size(); ++i) {
                    Point_3 &p = region_facet[i];

                    q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                    rotate_point(angle, axis, q);
                    p = Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z());
                }
            }

            void rotate_point(const FT angle, const Vector_3 &axis, Point_3 &p) const {

				const double tmp_angle = CGAL::to_double(angle);

				const FT c = static_cast<FT>(std::cos(tmp_angle));
				const FT s = static_cast<FT>(std::sin(tmp_angle));

				const FT C = FT(1) - c;

				const FT x = axis.x();
				const FT y = axis.y();
				const FT z = axis.z();

				p = Point_3((x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
					  		(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
					  		(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
			}

            void triangulate_region_facets(const Region_facets &region_facets, CDT &cdt) const {

				Vertex_handles vhs;
                insert_points(region_facets, cdt, vhs);
                
                Final_constraints final_vhs;
                update_constraints(region_facets, vhs, final_vhs);

                insert_constraints(final_vhs, cdt);
            }

            void insert_points(const Region_facets &region_facets, CDT &cdt, Vertex_handles &vhs) const {
                
                cdt.clear();
                vhs.clear();

                vhs.resize(region_facets.size());
				for (size_t i = 0; i < region_facets.size(); ++i) {
					const Region_facet &region_facet = region_facets[i];

					vhs[i].resize(region_facet.size());
					for (size_t j = 0; j < region_facet.size(); ++j) {
                        const Point_3 &p = region_facet[j];

						vhs[i][j] = cdt.insert(Point_2(p.x(), p.y()));
						vhs[i][j]->info().height = p.z();
					}
				}
            }

            void update_constraints(const Region_facets &region_facets, const Vertex_handles &vhs, Final_constraints &final_vhs) const {

                for (size_t i = 0; i < region_facets.size(); ++i) {
                    // new_vhs[i].resize(region_facets[i].size());

                    for (size_t j = 0; j < region_facets[i].size(); ++j) {
                        const size_t jp = (j + 1) % region_facets[i].size();

                        if (is_boundary_edge(region_facets[i][j], region_facets[i][jp], i, region_facets)) {
                            
                            // new_vhs[i][j]  = vhs[i][j];
                            // new_vhs[i][jp] = vhs[i][jp];

                            const Final_constraint final_constraint = std::make_pair(vhs[i][j], vhs[i][jp]);
                            final_vhs.push_back(final_constraint);
                            
                            // std::cout << "init constr: " << vhs[i][j]->point() << "; " << vhs[i][jp]->point() << std::endl;
                        }
                    }
                }
            }

            bool is_boundary_edge(const Point_3 &p1, const Point_3 &p2, const size_t facet_index, const Region_facets &region_facets) const {

                for (size_t i = 0; i < region_facets.size(); ++i) {
                    if (i == facet_index) continue;

                    for (size_t j = 0; j < region_facets[i].size(); ++j) {
                        const size_t jp = (j + 1) % region_facets[i].size();

                        if (are_equal_edges(p1, p2, region_facets[i][j], region_facets[i][jp])) return false;
                    }
                }
                return true;
            }

            bool are_equal_edges(const Point_3 &p1, const Point_3 &p2, const Point_3 &q1, const Point_3 &q2) const {
                return (are_equal_points(p1, q1) && are_equal_points(p2, q2)) || (are_equal_points(p1, q2) && are_equal_points(p2, q1));
            }

            template<class Point>
            bool are_equal_points(const Point &p, const Point &q) const {

                const FT eps = m_tolerance;
                return (CGAL::abs(p.x() - q.x()) < eps) && (CGAL::abs(p.y() - q.y()) < eps) && (CGAL::abs(p.z() - q.z()) < eps);
            }

            void insert_constraints(const Final_constraints &final_vhs, CDT &cdt) const {
                
                for (size_t i = 0; i < final_vhs.size(); ++i) {
                    const Final_constraint &final_constraint = final_vhs[i];
                    
                    if (final_constraint.first != final_constraint.second)
                        cdt.insert_constraint(final_constraint.first, final_constraint.second);

                    // std::cout << "insert constr: " << final_constraint.first->point() << "; " << final_constraint.second->point() << std::endl;
                }
            }

            void get_back_region_facets(const CDT &cdt, Region_facets &region_facets) const {

                if (m_use_all_facets) {
                
                    get_all_triangles(cdt, region_facets);

                } else get_only_boundary(cdt, region_facets);
            }

            void get_all_triangles(const CDT &cdt, Region_facets &region_facets) const {
                
                if (cdt.number_of_faces() < 2) {
                    if (cdt.number_of_faces() == 0) return;

                    /*
                    const Faces_iterator &fh = cdt.finite_faces_begin();

                    const Point_2 &p1 = fh->vertex(0)->point();
                    const Point_2 &p2 = fh->vertex(1)->point();
                    const Point_2 &p3 = fh->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (triangle.area() < m_area_tolerance) return; */
                }

                // std::cout << "cdt: " << cdt.number_of_faces() << std::endl;

                region_facets.clear();

                size_t i = 0; Region_facet region_facet(3);
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++i) {

                    const Vertex_handle &vh1 = fit->vertex(0);
                    const Vertex_handle &vh2 = fit->vertex(1);
                    const Vertex_handle &vh3 = fit->vertex(2);

                    const Point_2 &p1 = vh1->point();
                    const Point_2 &p2 = vh2->point();
                    const Point_2 &p3 = vh3->point();

                    // std::cout << p1 << "; " << p2 << "; " << p3 << std::endl;

                    if (vh1->info().height == m_default_height || 
                        vh2->info().height == m_default_height ||
                        vh3->info().height == m_default_height) continue;

                    region_facet[0] = Point_3(p1.x(), p1.y(), vh1->info().height);
                    region_facet[1] = Point_3(p2.x(), p2.y(), vh2->info().height);
                    region_facet[2] = Point_3(p3.x(), p3.y(), vh3->info().height);

                    region_facets.push_back(region_facet);
                }
            }

            void get_only_boundary(const CDT &cdt, Region_facets &region_facets) const {

                // std::cout << "cdt: " << cdt.number_of_faces() << std::endl;
                if (cdt.number_of_faces() == 0) return;

                Face_handle fh;
                Region_facet region_facet;

                bool success = find_first_face_handle(cdt, fh);
                if (!success) return;

                success = traverse_cdt(fh, cdt, region_facet);
                if (!success) return;

                if (region_facet.size() < 3) return;

                /*
                for (size_t i = 0; i < region_facet.size(); ++i)
                    std::cout << region_facet[i] << std::endl;
                std::cout << std::endl; */

                region_facets.clear();
                region_facets.push_back(region_facet);
            }

            bool find_first_face_handle(const CDT &cdt, Face_handle &fh) const {

                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    fh = static_cast<Face_handle>(fit);

                    const Vertex_handle &vh1 = fh->vertex(0);
                    const Vertex_handle &vh2 = fh->vertex(1);
                    const Vertex_handle &vh3 = fh->vertex(2);

                    const Point_2 &p1 = vh1->point();
                    const Point_2 &p2 = vh2->point();
                    const Point_2 &p3 = vh3->point();

                    // std::cout << p1 << "; " << p2 << "; " << p3 << std::endl;
                    for (size_t i = 0; i < 3; ++i) {
                        
                        const Edge edge = std::make_pair(fh, i);
                        if (cdt.is_constrained(edge)) {
                         
                            // std::cout << "constr: " << edge.first->vertex((edge.second + 1) %3)->point() << "; " << edge.first->vertex((edge.second + 2) %3)->point() << std::endl;
                            return true;
                        }
                    }
                }
                return false;
            }

            bool traverse_cdt(const Face_handle &fh, const CDT &cdt, Region_facet &region_facet) const {
                
                Edge edge;
                region_facet.clear();

                const bool success = find_first_edge(fh, cdt, edge);
                if (!success) return false;

                CGAL_precondition(edge.second >= 0 && edge.second <= 2);

                Vertex_handle  vh = edge.first->vertex((edge.second + 2) % 3);
                Vertex_handle end = vh;
                
                if (vh->info().height == m_default_height) return false;

                const Point_2 &p = vh->point();
                region_facet.push_back(Point_3(p.x(), p.y(), vh->info().height));

                // std::cout << "first p: " << p << std::endl;

                size_t num_iters = 0; 
                do {
                    get_next_vertex_handle(cdt, vh, edge);
                    const Point_2 &q = vh->point();

                    if (vh->info().height == m_default_height) return false;
                    if (vh == end) break;

                    // std::cout << "next q: " << q << std::endl;
                    region_facet.push_back(Point_3(q.x(), q.y(), vh->info().height));
                    
                    if (num_iters == m_max_num_iters) return false;
                    ++num_iters;

                } while (vh != end);

                return is_valid_traversal(region_facet);
            }

            bool is_valid_traversal(const Region_facet &region_facet) const {

                if (region_facet.size() < 3) return false;

                for (size_t i = 0; i < region_facet.size(); ++i) {
                    const Point_3 &p = region_facet[i];
                    
                    // std::cout << "p: " << p << std::endl;

                    for (size_t j = 0; j < region_facet.size(); ++j) {
                        const Point_3 &q = region_facet[j];

                        if (i == j) continue;

                        // std::cout << "q: " << q << std::endl;

                        if (are_equal_points(p, q)) return false;
                    }
                }
                // std::cout << std::endl;
                return true;
            }

            bool find_first_edge(const Face_handle &fh, const CDT &cdt, Edge &edge) const {

                for (int i = 0; i < 3; ++i) {
                    
                    edge = std::make_pair(fh, i);
                    if (cdt.is_constrained(edge)) {
                     
                        // std::cout << "first edge: " << edge.first->vertex((edge.second + 1) %3)->point() << "; " << edge.first->vertex((edge.second + 2) %3)->point() << std::endl;
                        return true;
                    }
                }
                return false;
            }

            void get_next_vertex_handle(const CDT &cdt, Vertex_handle &vh, Edge &edge) const {

                const int index = edge.first->index(vh);
                Edge next = std::make_pair(edge.first, (index + 2) % 3);

                while (!cdt.is_constrained(next)) {

                    const Face_handle fh = next.first->neighbor(next.second);
                    const Vertex_handle tmp = next.first->vertex((next.second + 1) % 3);

                    // std::cout << "neigh: " << fh->vertex(0)->point() << "; " << fh->vertex(1)->point() << "; " << fh->vertex(2)->point() << std::endl;
                    
                    const size_t tmp_index = fh->index(tmp);
                    next = std::make_pair(fh, (tmp_index + 2) % 3);
                }

                vh   = next.first->vertex((next.second + 2) % 3);
                edge = next;
            }

            Color generate_random_color() const {

				const int r = rand() % 255;
				const int g = rand() % 255;
				const int b = rand() % 255;

				return Color(r, g, b);
			}
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_MERGER_H