#ifndef CGAL_LEVEL_OF_DETAIL_FACET_REGIONS_MERGER_3_H
#define CGAL_LEVEL_OF_DETAIL_FACET_REGIONS_MERGER_3_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding>
		class Level_of_detail_facet_regions_merger_3 {

		public:
			using Kernel   = InputKernel;
            using Building = InputBuilding;

			using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Vector_3   = typename Kernel::Vector_3;

            using CDT = typename Building::CDT;

            using Region_facet   = std::vector<Point_3>;
            using Region_facets  = std::vector<Region_facet>;
            using Output_regions = std::vector<Region_facets>;

            using Clean_facets = Region_facets;

            typename Kernel::Compute_squared_length_3         squared_length_3;
            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Vertex_handle  = typename CDT::Vertex_handle;
            using Vertex_handles = std::vector< std::vector<Vertex_handle> >;

            using Face_handle    = typename CDT::Face_handle;
            using Faces_iterator = typename CDT::Finite_faces_iterator;
            
            using Edge = typename CDT::Edge;

            using Final_constraint  = std::pair<Vertex_handle, Vertex_handle>;
            using Final_constraints = std::vector<Final_constraint>;

            using Polygon_2 = CGAL::Polygon_2<Kernel>;

			Level_of_detail_facet_regions_merger_3(Output_regions &output_regions) :
            m_output_regions(output_regions),
            m_tolerance(FT(1) / FT(100000)),
            m_default_height(-FT(100000000000000)),
            m_use_all_facets(false),
            m_max_num_iters(100)
            { }

			void merge(Clean_facets &clean_facets) const {

                if (m_output_regions.size() == 0) return;
                merge_facets(m_output_regions, clean_facets);
            }

            void use_all_facets(const bool new_state) {
                m_use_all_facets = new_state;
            }

		private:
            Output_regions &m_output_regions;
            
            const FT m_tolerance;
            const FT m_default_height;
            
            bool m_use_all_facets;
            const size_t m_max_num_iters;

            void merge_facets(Output_regions &output_regions, Clean_facets &clean_facets) const {

                clean_facets.clear();
                const size_t out_size = output_regions.size();

                for (size_t i = 0; i < output_regions.size(); ++i)
                    merge_region_facets(output_regions[i], clean_facets);
            }

            bool merge_region_facets(Region_facets &region_facets, Clean_facets &clean_facets) const {

                Vector_3 source_normal;
                bool success = compute_source_normal(region_facets, source_normal);

                if (!success) {
                 
                    for (size_t i = 0; i < region_facets.size(); ++i)
                        clean_facets.push_back(region_facets[i]);
                    
                    return;
                }

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                if (source_normal == -target_normal) source_normal = target_normal;

                FT angle; Vector_3 axis;
                success = compute_angle_and_axis(source_normal, target_normal, angle, axis);

                if (!success) return;

                Point_3 b;
                compute_barycentre(region_facets, b);

                const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);
                
                if (angle_deg != FT(0) && angle_deg != FT(180))
                    rotate_region_facets(region_facets, angle, axis, b);

                CDT cdt;
                triangulate_region_facets(region_facets, cdt);
                
                if (cdt.number_of_faces() != 0) get_back_region_facets(cdt, region_facets);
                fix_orientation(region_facets);

                if (angle_deg != FT(0) && angle_deg != FT(180))
                    rotate_region_facets(region_facets, -angle, axis, b, true);
                
                for (size_t i = 0; i < region_facets.size(); ++i)
                    clean_facets.push_back(region_facets[i]);
            }

            bool compute_source_normal(const Region_facets &region_facets, Vector_3 &normal) const {
                
                if (region_facets.size() == 0) 
                    return false;

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

                if (are_equal_points(normal, zero)) 
                    return false;

                normalize(normal);
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
                    if (!are_equal_points(normal, zero)) return true;
                }
                return false;

                std::cout << "error: normal not found" << std::endl;
                exit(1);
            }

            template<class Point>
            bool are_equal_points(const Point &p, const Point &q) const {

                const FT eps = m_tolerance;
                return (CGAL::abs(p.x() - q.x()) < eps) && (CGAL::abs(p.y() - q.y()) < eps) && (CGAL::abs(p.z() - q.z()) < eps);
            }

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void compute_target_normal(Vector_3 &normal) const {
                normal = Vector_3(FT(0), FT(0), FT(1));
            }

            bool compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {

				const auto  cross = cross_product_3(m, n);
				const   FT length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const   FT    dot = dot_product_3(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
				const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);

				if (angle_deg == FT(0) || angle_deg == FT(180)) 
					return true;

				if (length == FT(0)) {
                 
                    std::cout << "error weight quality estimator: length = 0" << std::endl;
                    exit(0);
                 
                    return false;
                }
                
				CGAL_precondition(length != FT(0));
				axis = cross / length;

                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle > half_pi) {
                    
                    angle = static_cast<FT>(CGAL_PI) - angle;
                    axis = -axis;
                }

				return true;
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

            void rotate_region_facets(Region_facets &region_facets, const FT angle, const Vector_3 &axis, const Point_3 &b, const bool show = false) const {

                for (size_t i = 0; i < region_facets.size(); ++i)
                    rotate_region_facet(region_facets[i], angle, axis, b);
            }

            void rotate_region_facet(Region_facet &region_facet, const FT angle, const Vector_3 &axis, const Point_3 &b) const {

                if (angle != FT(0)) 
                    rotate_facet(b, angle, axis, region_facet);
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

                    for (size_t j = 0; j < region_facets[i].size(); ++j) {
                        const size_t jp = (j + 1) % region_facets[i].size();

                        if (is_boundary_edge(region_facets[i][j], region_facets[i][jp], i, region_facets)) {

                            const Final_constraint final_constraint = std::make_pair(vhs[i][j], vhs[i][jp]);
                            final_vhs.push_back(final_constraint);
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

            void insert_constraints(const Final_constraints &final_vhs, CDT &cdt) const {
                
                for (size_t i = 0; i < final_vhs.size(); ++i) {
                    const Final_constraint &final_constraint = final_vhs[i];
                    
                    if (final_constraint.first != final_constraint.second)
                        cdt.insert_constraint(final_constraint.first, final_constraint.second);
                }
            }

            void get_back_region_facets(const CDT &cdt, Region_facets &region_facets) const {

                if (m_use_all_facets) get_all_triangles(cdt, region_facets);
                else get_only_boundary(cdt, region_facets);
            }

            void get_all_triangles(const CDT &cdt, Region_facets &region_facets) const {
                
                if (cdt.number_of_faces() < 2)
                    if (cdt.number_of_faces() == 0) return;

                region_facets.clear();

                size_t i = 0; Region_facet region_facet(3);
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++i) {

                    const Vertex_handle &vh1 = fit->vertex(0);
                    const Vertex_handle &vh2 = fit->vertex(1);
                    const Vertex_handle &vh3 = fit->vertex(2);

                    const Point_2 &p1 = vh1->point();
                    const Point_2 &p2 = vh2->point();
                    const Point_2 &p3 = vh3->point();

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

                if (cdt.number_of_faces() == 0) return;

                Face_handle fh;
                Region_facet region_facet;

                bool success = find_first_face_handle(cdt, fh);
                if (!success) return;

                success = traverse_cdt(fh, cdt, region_facet);
                if (!success) return;

                if (region_facet.size() < 3) return;

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

                    for (size_t i = 0; i < 3; ++i) {
                        
                        const Edge edge = std::make_pair(fh, i);
                        if (cdt.is_constrained(edge)) return true;
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

                size_t num_iters = 0; 
                do {
                    get_next_vertex_handle(cdt, vh, edge);
                    const Point_2 &q = vh->point();

                    if (vh->info().height == m_default_height) return false;
                    if (vh == end) break;

                    region_facet.push_back(Point_3(q.x(), q.y(), vh->info().height));
                    
                    if (num_iters == m_max_num_iters) return false;
                    ++num_iters;

                } while (vh != end);

                return is_valid_traversal(region_facet);
            }

            bool find_first_edge(const Face_handle &fh, const CDT &cdt, Edge &edge) const {

                for (int i = 0; i < 3; ++i) {
                    
                    edge = std::make_pair(fh, i);
                    if (cdt.is_constrained(edge)) return true;
                }
                return false;
            }

            void get_next_vertex_handle(const CDT &cdt, Vertex_handle &vh, Edge &edge) const {

                const int index = edge.first->index(vh);
                Edge next = std::make_pair(edge.first, (index + 2) % 3);

                while (!cdt.is_constrained(next)) {

                    const Face_handle fh = next.first->neighbor(next.second);
                    const Vertex_handle tmp = next.first->vertex((next.second + 1) % 3);
                    
                    const size_t tmp_index = fh->index(tmp);
                    next = std::make_pair(fh, (tmp_index + 2) % 3);
                }

                vh   = next.first->vertex((next.second + 2) % 3);
                edge = next;
            }

            bool is_valid_traversal(const Region_facet &region_facet) const {

                if (region_facet.size() < 3) return false;

                for (size_t i = 0; i < region_facet.size(); ++i) {
                    const Point_3 &p = region_facet[i];

                    for (size_t j = 0; j < region_facet.size(); ++j) {
                        const Point_3 &q = region_facet[j];

                        if (i == j) continue;
                        if (are_equal_points(p, q)) return false;
                    }
                }
                return true;
            }

            void fix_orientation(Region_facets &region_facets) const {

                Polygon_2 polygon;
                for (size_t i = 0; i < region_facets.size(); ++i) {
                    Region_facet &region_facet = region_facets[i];

                    polygon.clear();
                    for (size_t j = 0; j < region_facet.size(); ++j) {
                        
                        const Point_3 &p = region_facet[j];
                        polygon.push_back(Point_2(p.x(), p.y()));
                    }

                    if (polygon.is_clockwise_oriented()) 
                        std::reverse(region_facet.begin(), region_facet.end());
                }
            }
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_FACET_REGIONS_MERGER_3_H