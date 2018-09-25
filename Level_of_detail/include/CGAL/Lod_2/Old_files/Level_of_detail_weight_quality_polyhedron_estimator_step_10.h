#ifndef CGAL_LEVEL_OF_DETAIL_WEIGHT_QUALITY_POLYHEDRON_ESTIMATOR_STEP_10_H
#define CGAL_LEVEL_OF_DETAIL_WEIGHT_QUALITY_POLYHEDRON_ESTIMATOR_STEP_10_H

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_weight_quality_polyhedron_estimator_step_10 {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Plane_3    = typename Kernel::Plane_3;
            using Vector_3   = typename Kernel::Vector_3;
            using Triangle_2 = typename Kernel::Triangle_2;
            
            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Index   = typename Building::Index;
            using Indices = typename Building::Indices;

            using Polyhedron_facet    = typename Polyhedron::Facet;
            using Polyhedron_facets   = typename Polyhedron::Facets;
            using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Graphcut_facet  = typename Building::Graphcut_facet;
            using Graphcut_facets = typename Building::Graphcut_facets;

            using Graphcut_facet_data      = typename Graphcut_facet::Data;
            using Graphcut_facet_data_pair = typename Graphcut_facet::Data_pair;

            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;

            using Wall_boundary = typename Wall::Wall_boundary;

            using Points_2 = std::vector<Point_2>;
            using Points_3 = std::vector<Point_3>;

            using Polygon = CGAL::Polygon_2<Kernel>;

            typename Kernel::Compute_squared_length_3 squared_length_3;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;
            typename Kernel::Compute_squared_distance_3 squared_distance_3;
            
            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Coordinates = std::vector<FT>;

            using Mean_value = CGAL::Barycentric_coordinates::Mean_value_2<Kernel>;
            using Mean_value_coordinates = CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel>;

            using Label     = int;
			using Label_map = typename Input:: template Property_map<Label>;

            using Log = CGAL::LOD::Mylog;

            using Vb  = CGAL::Alpha_shape_vertex_base_2<Kernel>;
            using Fb  = CGAL::Alpha_shape_face_base_2<Kernel>;
            using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
            
            using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
            using Alpha_shape_2   = CGAL::Alpha_shape_2<Triangulation_2>;

            using Faces_iterator = typename Triangulation_2::Finite_faces_iterator;
            using Face_handle    = typename Triangulation_2::Face_handle;

            Level_of_detail_weight_quality_polyhedron_estimator_step_10(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height),
            m_tolerance(FT(1) / FT(100000)),
            m_edge_length_tolerance(FT(1) / FT(1000)),
            m_distance_threshold(FT(2)),
            m_big_value(FT(100000000000000)),
            m_angle_threshold(FT(10)),
            m_extra_bc(FT(1) / FT(1000)),
            m_bc_tolerance_top(FT(1) + m_extra_bc),
            m_bc_tolerance_bottom(FT(0) - m_extra_bc),
            m_default_weight(FT(0)),
            m_default_quality(FT(0)),
            m_alpha(-FT(1)) { 

                boost::tie(m_labels, boost::tuples::ignore) = m_input.template property_map<Label>("label");
            }

			void set_alpha(const FT new_value) {

				CGAL_precondition(new_value > FT(0));
				m_alpha = new_value;
			}

            void estimate() {

                if (m_buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid && building.interior_indices.size() != 0)
                        process_building(building);
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;

            const FT m_ground_height;
            const FT m_tolerance;
            const FT m_edge_length_tolerance;
            const FT m_distance_threshold;
            const FT m_big_value;
            const FT m_angle_threshold;

            const FT m_extra_bc;
            const FT m_bc_tolerance_top;
            const FT m_bc_tolerance_bottom;

            Label_map m_labels;

            const FT m_default_weight;
            const FT m_default_quality;

            FT m_alpha;

            void process_building(Building &building) {

                Graphcut_facets &gc_facets = building.graphcut_facets;
                create_graphcut_facets(building, gc_facets);
            }

            void create_graphcut_facets(const Building &building, Graphcut_facets &gc_facets) const {
                
                gc_facets.clear();
                
                const Polyhedrons &polyhedrons = building.polyhedrons;

                Indices interior_indices;
                create_interior_indices(building, interior_indices);

                for (int i = 0; i < polyhedrons.size(); ++i)
                    add_graphcut_facets(interior_indices, polyhedrons, i, gc_facets);
            }

            void create_interior_indices(const Building &building, Indices &interior_indices) const {

                const Label facade = 1;
                const Label roof   = 2;

                interior_indices.clear();

                for (size_t i = 0; i < m_input.number_of_points(); ++i) {
                    const Point_3 &p = m_input.point(i);
        
                    if (m_labels[i] == facade || m_labels[i] == roof) {

                        if (belongs_to_building(p, building))
                            interior_indices.push_back(i);
                    }
                }
            }

            bool belongs_to_building(const Point_3 &p, const Building &building) const {

                Point_3 minp, maxp;
                compute_3d_bounding_box(building, minp, maxp);

                return belongs_to_bbox(p, minp, maxp);
            }

            void compute_3d_bounding_box(const Building &building, Point_3 &minp, Point_3 &maxp) const {

                const Walls &walls = building.walls;

                FT minx =  m_big_value, miny =  m_big_value, minz =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value, maxz = -m_big_value;

                for (size_t i = 0; i < walls.size(); ++i) {
                    const Wall_boundary &boundary = walls[i].boundary;
                    
                    for (size_t j = 0; j < boundary.size(); ++j) {
                        const Point_3 &p = boundary[j];

                        minx = CGAL::min(minx, p.x());
                        miny = CGAL::min(miny, p.y());

                        maxx = CGAL::max(maxx, p.x());
                        maxy = CGAL::max(maxy, p.y());
                    }
                }

                const FT extra_distance = m_distance_threshold;

                minx -= extra_distance;
                maxx += extra_distance;

                miny -= extra_distance;
                maxy += extra_distance;

                minz = m_ground_height;
                maxz = m_big_value;

                minp = Point_3(minx, miny, minz);
                maxp = Point_3(maxx, maxy, maxz);
            }

            bool belongs_to_bbox(const Point_3 &q, const Point_3 &minp, const Point_3 &maxp) const {
                return (q.x() > minp.x() && q.x() < maxp.x()) && (q.y() > minp.y() && q.y() < maxp.y()) && (q.z() > minp.z() && q.z() < maxp.z());
            }

            void add_graphcut_facets(const Indices &interior_indices, const Polyhedrons &polyhedrons, const int poly_index, Graphcut_facets &gc_facets) const {

                for (int i = 0; i < polyhedrons[poly_index].facets.size(); ++i)
                    process_polyhedron_facet(interior_indices, polyhedrons, poly_index, i, gc_facets);
            }

            void process_polyhedron_facet(const Indices &interior_indices, const Polyhedrons &polyhedrons, const int poly_index, const int facet_index, Graphcut_facets &gc_facets) const {

                if (was_already_added(poly_index, facet_index, gc_facets))
                    return;

                // std::cout << poly_index << " " << facet_index << std::endl;

                // if (poly_index  == 2 && facet_index == 1) {
                // if (poly_index == 12 && facet_index == 5) {

                const Polyhedron &polyhedron = polyhedrons[poly_index];
                    
                const Polyhedron_facets   &facets   = polyhedron.facets;
                const Polyhedron_vertices &vertices = polyhedron.vertices;

                Graphcut_facet gc_facet;
                gc_facet.data = std::make_pair(poly_index, facet_index);

                find_neighbours(polyhedrons, facets, vertices, poly_index, facet_index, gc_facet.neighbours);
                    
                compute_weight(facets[facet_index], vertices, gc_facet);
                compute_quality(interior_indices, facets[facet_index], vertices, gc_facet);

                if (gc_facet.quality > FT(2) / FT(5)) gc_facet.is_valid = true;
                else gc_facet.is_valid = false;

                gc_facets.push_back(gc_facet);
                    
                // exit(0);
                // }
            }

            bool was_already_added(const int poly_index, const int facet_index, const Graphcut_facets &gc_facets) const {

                for (size_t i = 0; i < gc_facets.size(); ++i) {
                    const Graphcut_facet_data_pair &data_pair = gc_facets[i].neighbours;

                    const Graphcut_facet_data &neigh_1 = data_pair.first;
                    const Graphcut_facet_data &neigh_2 = data_pair.second;

                    if (neigh_1.first == poly_index && neigh_1.second == facet_index) return true;
                    if (neigh_2.first == poly_index && neigh_2.second == facet_index) return true;
                }
                return false;
            }

            void compute_weight(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Graphcut_facet &gc_facet) const {

                const FT weight = compute_facet_area(facet, vertices);
                gc_facet.weight = get_weight(weight);
            }

            FT compute_facet_area(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices) const {

                Vector_3 source_normal;
                bool success = compute_source_normal(facet, vertices, source_normal);

                if (!success) return get_default_weight();

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                if (source_normal == -target_normal) source_normal = target_normal;

                FT angle; Vector_3 axis;
                success = compute_angle_and_axis(source_normal, target_normal, angle, axis);

                if (!success) return get_default_weight();
                
                Point_3 b;
                compute_3d_polygon_barycentre(facet, vertices, b);
                    
                Polygon polygon;
                success = create_polygon(facet, vertices, b, angle, axis, polygon);

                if (!success) return get_default_weight();
                return CGAL::abs(polygon.area());
            }

            FT get_default_weight() const {
                return m_default_weight;
            }

            FT get_weight(const FT weight) const {
                return weight;
            }

            void compute_quality(const Indices &interior_indices, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Graphcut_facet &gc_facet) const {

                const Graphcut_facet_data &neigh_1 = gc_facet.neighbours.first;
                const Graphcut_facet_data &neigh_2 = gc_facet.neighbours.second;

                if (neigh_1.first < 0 || neigh_1.second < 0) {
                    
                    gc_facet.quality = get_default_quality();
                    return;
                }

                if (neigh_2.first < 0 || neigh_2.second < 0) {
                
                    gc_facet.quality = get_default_quality();
                    return;    
                }

                if (is_vertical_facet(vertices, facet)) {

                    gc_facet.quality = FT(9) / FT(10);
                    return;   
                }

                return;

                if (is_ground_facet(vertices, facet)) {

                    gc_facet.quality = FT(9) / FT(10);
                    return;
                }

                Indices closest;
                find_closest_points(interior_indices, facet, vertices, closest);

                if (closest.size() == 0) {
                    
                    gc_facet.quality = get_default_quality();
                    return;
                }

                Indices internal;
                const FT initial_quality = find_internal_points(closest, facet, vertices, internal);

                if (internal.size() == 0) {
                    
                    gc_facet.quality = get_default_quality();
                    return;
                }

                const FT extra_quality = compute_extra_quality(closest, internal);
                gc_facet.quality       = compute_final_quality(initial_quality, extra_quality);

                // std::cout << internal.size() << " : " << closest.size() << " : " << initial_quality << " : " << extra_quality << " : " << gc_facet.quality << std::endl;
            }

            FT get_default_quality() const {
                return m_default_quality;
            }

            bool is_vertical_facet(const Polyhedron_vertices &vertices, const Polyhedron_facet &facet) const {

				Vector_3 facet_normal;
				set_facet_normal(vertices, facet, facet_normal);

				Vector_3 ground_normal;
				set_ground_normal(ground_normal);

                const FT angle      = compute_angle(facet_normal, ground_normal);
                const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle));

                if (angle_diff < m_angle_threshold) return true;
                return false;
            }

            void set_facet_normal(const Polyhedron_vertices &vertices, const Polyhedron_facet &facet, Vector_3 &facet_normal) const {

                CGAL_precondition(facet.indices.size() >= 3);
                const Point_3 &p1 = vertices[facet.indices[0]];
                const Point_3 &p2 = vertices[facet.indices[1]];
                const Point_3 &p3 = vertices[facet.indices[2]];

                const Vector_3 v1 = Vector_3(p1, p2);
                const Vector_3 v2 = Vector_3(p1, p3);

                facet_normal = cross_product_3(v1, v2);
                normalize(facet_normal);
			}

			void set_ground_normal(Vector_3 &ground_normal) const {
				ground_normal = Vector_3(FT(0), FT(0), FT(1));
			}

            FT compute_angle(const Vector_3 &m, const Vector_3 &n) const {
				
				const auto cross = cross_product_3(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const FT dot     = dot_product_3(m, n);

				FT angle_rad = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                
                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle_rad > half_pi) angle_rad = static_cast<FT>(CGAL_PI) - angle_rad;

				const FT angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
                return angle_deg;
			}

            bool is_ground_facet(const Polyhedron_vertices &vertices, const Polyhedron_facet &facet) const {

                FT average_height = FT(0);
                for (size_t i = 0; i < facet.indices.size(); ++i) {
                    
                    const Point_3 &p = vertices[facet.indices[i]];
                    average_height += p.z();
                }
                average_height /= static_cast<FT>(facet.indices.size());

                if (CGAL::abs(average_height - m_ground_height) < m_distance_threshold / FT(2)) return true;
                return false;
            }

            void find_closest_points(const Indices &interior_indices, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Indices &closest) const {

                Plane_3 plane;
                const bool success = create_plane(facet, vertices, plane);

                if (!success) return;
                find_closest_points_to_plane(interior_indices, plane, closest);
            }

            bool create_plane(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Plane_3 &plane) const {

                if (facet.indices.size() < 3) return false;
                
                Points_3 points;
                get_base_points(facet, vertices, points);
                
                if (points.size() < 3) return false;

                const Point_3 &p1 = points[0];
                const Point_3 &p2 = points[1];
                const Point_3 &p3 = points[2];

                plane = Plane_3(p1, p2, p3);
                return true;
            }

            void get_base_points(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Points_3 &points) const {

                points.clear();
                const size_t num_indices = facet.indices.size();

                for (size_t i = 0; i < num_indices; ++i) {
                    const size_t ip = (i + 1) % num_indices;

                    const Point_3 &p = vertices[facet.indices[i]];
                    const Point_3 &q = vertices[facet.indices[ip]];

                    const FT edge_length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_3(p, q))));

                    if (edge_length < m_edge_length_tolerance) continue;
                    else points.push_back(p);
                }
            }

            void find_closest_points_to_plane(const Indices &interior_indices, const Plane_3 &plane, Indices &closest) const {

                closest.clear();
                
                if (interior_indices.size() == 0) 
                    return;

                for (size_t i = 0; i < interior_indices.size(); ++i) {
                    
                    const Point_3 &p = m_input.point(interior_indices[i]);
                    const Point_3  q = plane.projection(p);

                    const FT squared_distance = squared_distance_3(p, q);
                    const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance)));

                    if (distance < m_distance_threshold)
                        closest.push_back(interior_indices[i]);
                }
            }

            FT find_internal_points(const Indices &closest, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Indices &internal) const {

                internal.clear();
                if (closest.size() == 0) return get_default_quality();

                Vector_3 source_normal;
                bool success = compute_source_normal(facet, vertices, source_normal);

                if (!success) return get_default_quality();

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                if (source_normal == -target_normal) source_normal = target_normal;

                FT angle; Vector_3 axis;
                success = compute_angle_and_axis(source_normal, target_normal, angle, axis);

                if (!success) return get_default_quality();

                Point_3 b;
                compute_3d_polygon_barycentre(facet, vertices, b);
                    
                Polygon polygon;
                success = create_polygon(facet, vertices, b, angle, axis, polygon);

                if (!success) return get_default_quality();

                Points_2 queries;
                create_queries(closest, b, angle, axis, queries);

                // debug

                /*

                Log logger;
                logger.export_polygon("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_polygon", facet, vertices);
                logger.export_polygon("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_polygon_rotated", polygon);
                
                Points_3 closest_points(closest.size());
                for (size_t i = 0; i < closest.size(); ++i)
                    closest_points[i] = m_input.point(closest[i]);

                logger.export_3d_points("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_points", closest_points);
                logger.export_points("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_points_rotated", queries);

                */

                //

                Indices tmp;
                find_points_inside_polygon(polygon, queries, closest, tmp, internal);
                
                if (internal.size() == 0) 
                    return get_default_quality();

                return estimate_initial_quality(polygon, queries, tmp);
            }

            bool compute_source_normal(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Vector_3 &normal) const {

                if (facet.indices.size() < 3) return false;
                
                Points_3 points;
                get_base_points(facet, vertices, points);
                
                if (points.size() < 3) return false;

                const Point_3 &p1 = points[0];
                const Point_3 &p2 = points[1];
                const Point_3 &p3 = points[2];

                const Vector_3 v1 = Vector_3(p2, p1);
                const Vector_3 v2 = Vector_3(p2, p3);

                normal = cross_product_3(v1, v2);
                
                const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
                if (are_equal_points(normal, zero)) return false;

                normalize(normal);
                return true;
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

            void compute_3d_polygon_barycentre(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Point_3 &b) const {
                
                CGAL_precondition(facet.indices.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < facet.indices.size(); ++i) {
                    const Point_3 &p = vertices[facet.indices[i]];

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(facet.indices.size());
                y /= static_cast<FT>(facet.indices.size());
                z /= static_cast<FT>(facet.indices.size());

                b = Point_3(x, y, z);
            }

            /*
            void compute_barycentre(const Indices &indices, Point_3 &b) const {
                
                CGAL_precondition(indices.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < indices.size(); ++i) {
                    const Point_3 &p = m_input.point(indices[i]);

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(indices.size());
                y /= static_cast<FT>(indices.size());
                z /= static_cast<FT>(indices.size());

                b = Point_3(x, y, z);
            } */

            bool create_polygon(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, const Point_3 &b, const FT angle, const Vector_3 &axis, Polygon &polygon) const {

                polygon.clear();

                for (size_t i = 0; i < facet.indices.size(); ++i) {
                    const Point_3 &p = vertices[facet.indices[i]];

				    const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);
                    if (angle_deg != FT(0) && angle_deg != FT(180)) {

                        Point_3 q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                        rotate_point(angle, axis, q);
                        polygon.push_back(Point_2(q.x() + b.x(), q.y() + b.y()));

                    } else polygon.push_back(Point_2(p.x(), p.y()));
                }

                if (!polygon.is_simple()) 
                    return false;

                if (polygon.is_clockwise_oriented()) 
                    polygon.reverse_orientation();

                return polygon.size() >= 3;
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

            void create_queries(const Indices &closest, const Point_3 &b, const FT angle, const Vector_3 &axis, Points_2 &queries) const {

                queries.clear();

                for (size_t i = 0; i < closest.size(); ++i) {
                    const Point_3 &p = m_input.point(closest[i]);

                    const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);
                    if (angle_deg != FT(0) && angle_deg != FT(180)) {

                        Point_3 q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                        rotate_point(angle, axis, q);
                        queries.push_back(Point_2(q.x() + b.x(), q.y() + b.y()));

                    } else queries.push_back(Point_2(p.x(), p.y()));
                }
            }

            void find_points_inside_polygon(const Polygon &polygon, const Points_2 &queries, const Indices &closest, Indices &tmp, Indices &internal) const {
                
                CGAL_precondition(queries.size() == closest.size());

                tmp.clear();
                internal.clear();

                Mean_value_coordinates mean_value_coordinates(polygon.vertices_begin(), polygon.vertices_end());

                for (size_t i = 0; i < queries.size(); ++i) {
                    const Point_2 &query = queries[i];

                    Coordinates coordinates;
                    mean_value_coordinates(query, std::back_inserter(coordinates));

                    if (is_inside_polygon(coordinates)) {
                     
                        tmp.push_back(i);
                        internal.push_back(closest[i]);
                    }
                }

                // std::cout << internal.size() << " : " << tmp.size() << std::endl;
            }

            bool is_inside_polygon(const Coordinates &coordinates) const {
                
                for (size_t i = 0 ; i < coordinates.size(); ++i)
                    if (coordinates[i] < m_bc_tolerance_bottom || coordinates[i] > m_bc_tolerance_top) return false;
                return true;
            }

            FT estimate_initial_quality(const Polygon &polygon, const Points_2 &queries, const Indices &tmp) const {
                return compute_alpha_shapes_based_quality(polygon, queries, tmp);
            }

            FT compute_alpha_shapes_based_quality(const Polygon &polygon, const Points_2 &queries, const Indices &tmp) const {

                const FT area1 = compute_polygon_area(polygon);
                const FT area2 = compute_alpha_shapes_area(queries, tmp);

                FT quality = -FT(1);

                if (area1 > area2) quality = area2 / area1;
                else quality = area1 / area2;

                CGAL_precondition(quality >= FT(0) && quality <= FT(1));

                const FT initial_quality = get_quality(quality);
                return initial_quality;
            }

            FT compute_polygon_area(const Polygon &polygon) const {

                // debug

                /*
                Log logger;
                logger.export_polygon("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "in_polygon", polygon);
                std::cout << "polygon area: " << polygon.area() << std::endl; */

                //

                return polygon.area();
            }

            FT compute_alpha_shapes_area(const Points_2 &queries, const Indices &tmp) const {

				CGAL_precondition(queries.size() != 0);
                CGAL_precondition(tmp.size()     != 0);

				Triangulation_2 triangulation;
                for (size_t i = 0; i < tmp.size(); ++i) {
                    
                    const Point_2 &p = queries[tmp[i]];
                    triangulation.insert(p);
                }

				CGAL_precondition(m_alpha > FT(0));
				Alpha_shape_2 in_cdt(triangulation, m_alpha, Alpha_shape_2::GENERAL);

                FT total_area = FT(0);
                for (Faces_iterator fit = in_cdt.finite_faces_begin(); fit != in_cdt.finite_faces_end(); ++fit) {
                    const Face_handle &fh = static_cast<Face_handle>(fit);

                    const typename Alpha_shape_2::Classification_type cl_type = in_cdt.classify(fh);
                    if (cl_type == Alpha_shape_2::INTERIOR) {

                        const Point_2 &p1 = fh->vertex(0)->point();
                        const Point_2 &p2 = fh->vertex(1)->point();
                        const Point_2 &p3 = fh->vertex(2)->point();

                        const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                        total_area += triangle.area();
                    }
                }

                // std::cout << "alpha shapes area: " << total_area << std::endl << std::endl;
                
                return total_area;

                // debug

                /*

                Points_2 in_points;
                for (size_t i = 0; i < tmp.size(); ++i) {
                    
                    const Point_2 &p = queries[tmp[i]];
                    in_points.push_back(p);
                }

                Log logger;
                logger.export_points("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "in_points", in_points);
                logger.save_cdt_ply(in_cdt, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "in_cdt");

                Points_2 out_points;
				for (auto ait = in_cdt.alpha_shape_vertices_begin(); ait != in_cdt.alpha_shape_vertices_end(); ++ait)
					out_points.push_back((*ait)->point());

                std::vector<Triangle_2> out_cdt;
                for (Faces_iterator fit = in_cdt.finite_faces_begin(); fit != in_cdt.finite_faces_end(); ++fit) {
                    const Face_handle &fh = static_cast<Face_handle>(fit);

                    const typename Alpha_shape_2::Classification_type cl_type = in_cdt.classify(fh);
                    if (cl_type == Alpha_shape_2::INTERIOR) {

                        const Point_2 &p1 = fh->vertex(0)->point();
                        const Point_2 &p2 = fh->vertex(1)->point();
                        const Point_2 &p3 = fh->vertex(2)->point();

                        const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                        out_cdt.push_back(triangle);
                    }
                }

                logger.export_points("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "out_points", out_points);
                logger.save_triangles(out_cdt, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "out_cdt");

                exit(0);

                // */
            }

            FT get_quality(const FT quality) const {
                
                CGAL_precondition(quality >= FT(0) && quality <= FT(1));
                return quality;
            }

            /*
            FT compute_bbox_based_quality(const Polygon &polygon, const Points_2 &queries, const Indices &tmp) const {

                Points_2 bbox1;
                compute_2d_polygon_bounding_box(polygon, bbox1);

                Points_2 bbox2;
                compute_2d_points_bounding_box(queries, tmp, bbox2);

                const FT area1 = compute_bbox_area(bbox1);
                const FT area2 = compute_bbox_area(bbox2);

                FT quality = -FT(1);

                if (area1 > area2) quality = area2 / area1;
                else quality = area1 / area2;

                CGAL_precondition(quality >= FT(0) && quality <= FT(1));

                const FT initial_quality = get_quality(quality);
                return initial_quality;
            }

            FT inverse(const FT value) const {
                
                CGAL_precondition(value >= FT(0) && value <= FT(1));
                return FT(1) - value;
            }

            void compute_2d_polygon_bounding_box(const Polygon &polygon, Points_2 &bbox) const {
                
                CGAL_precondition(polygon.size() > 2);
                
                FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;

                for (size_t i = 0; i < polygon.size(); ++i) {
                    const Point_2 &p = polygon.vertex(i);

                    minx = CGAL::min(minx, p.x());
                    miny = CGAL::min(miny, p.y());

                    maxx = CGAL::max(maxx, p.x());
                    maxy = CGAL::max(maxy, p.y());
                }

                bbox.clear();
                bbox.resize(4);

                bbox[0] = Point_2(minx, miny);
                bbox[1] = Point_2(maxx, miny);
                bbox[2] = Point_2(maxx, maxy);
                bbox[3] = Point_2(minx, maxy);
            }

            void compute_2d_points_bounding_box(const Points_2 &queries, const Indices &tmp, Points_2 &bbox) const {
                
                CGAL_precondition(tmp.size() > 0);
                
                FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;

                for (size_t i = 0; i < tmp.size(); ++i) {
                    const Point_2 &p = queries[tmp[i]];

                    minx = CGAL::min(minx, p.x());
                    miny = CGAL::min(miny, p.y());

                    maxx = CGAL::max(maxx, p.x());
                    maxy = CGAL::max(maxy, p.y());
                }

                bbox.clear();
                bbox.resize(4);

                bbox[0] = Point_2(minx, miny);
                bbox[1] = Point_2(maxx, miny);
                bbox[2] = Point_2(maxx, maxy);
                bbox[3] = Point_2(minx, maxy);
            }

            FT compute_bbox_area(const Points_2 &bbox) const {

                const FT width  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(bbox[0], bbox[1]))));
                const FT height = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(bbox[1], bbox[2]))));

                return width * height;
            }

            void compute_barycentre(const Polygon &polygon, Point_2 &b) const {

                CGAL_precondition(polygon.size() > 2);
                FT x = FT(0), y = FT(0);

                for (size_t i = 0; i < polygon.size(); ++i) {
                    const Point_2 &p = polygon.vertex(i);

                    x += p.x();
                    y += p.y();
                }

                x /= static_cast<FT>(polygon.size());
                y /= static_cast<FT>(polygon.size());

                b = Point_2(x, y);
            }

            void compute_barycentre(const Points_2 &queries, const Indices &tmp, Point_2 &b) const {

                CGAL_precondition(tmp.size() != 0);
                FT x = FT(0), y = FT(0);

                for (size_t i = 0; i < tmp.size(); ++i) {
                    const Point_2 &p = queries[tmp[i]];

                    x += p.x();
                    y += p.y();
                }

                x /= static_cast<FT>(tmp.size());
                y /= static_cast<FT>(tmp.size());

                b = Point_2(x, y);
            }

            void find_furthest_vertex(const Polygon &polygon, const Point_2 &q, Point_2 &b) const {

                CGAL_precondition(polygon.size() > 2);
                
                int index = -1; FT max_distance = -m_big_value;
                for (int i = 0; i < polygon.size(); ++i) {
                    
                    const Point_2 &p = polygon.vertex(i);
                    const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(p, q))));

                    if (distance > max_distance) {

                        index = i;
                        max_distance = distance;
                    }
                }

                CGAL_postcondition(index >= 0);
                b = polygon.vertex(index);
            } */

            FT compute_extra_quality(const Indices &closest, const Indices &internal) const {

                CGAL_precondition(internal.size() <= closest.size());

                const FT closest_size  = static_cast<FT>(closest.size());
                const FT internal_size = static_cast<FT>(internal.size());

                const FT quality = internal_size / closest_size;
                const FT extra_quality = get_quality(quality);

                return extra_quality;
            }

            FT compute_final_quality(const FT initial_quality, const FT extra_quality) const {

                CGAL_precondition(initial_quality >= FT(0) && initial_quality <= FT(1));
                CGAL_precondition(extra_quality   >= FT(0) && extra_quality   <= FT(1));

                const FT final_quality = initial_quality; // CGAL::min(initial_quality, extra_quality);
                return get_final_quality(final_quality);
            }

            FT get_final_quality(const FT final_quality) const {

                return const_function_value(final_quality);

                // return tanh_function_value(final_quality);
            }

            FT const_function_value(const FT value) const {
                return value;
            }

            FT tanh_function_value(const FT value) const {
                return static_cast<FT>(std::tanh(4.0 * CGAL::to_double(value)));
            }

            void find_neighbours(const Polyhedrons &polyhedrons, const Polyhedron_facets &facets, const Polyhedron_vertices &vertices, 
            const int poly_index, const int facet_index, Graphcut_facet_data_pair &neighbours) const {

                Graphcut_facet_data &neigh_1 = neighbours.first;
                Graphcut_facet_data &neigh_2 = neighbours.second;

                neigh_1.first  = poly_index;
                neigh_1.second = facet_index;

                find_neighbour(polyhedrons, poly_index, facets[facet_index], vertices, neigh_2);
            }

            void find_neighbour(const Polyhedrons &polyhedrons, const int poly_index, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Graphcut_facet_data &neigh) const {

                neigh.first  = -1;
                neigh.second = -1;

                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    if (i == poly_index) continue;

                    for (size_t j = 0; j < polyhedrons[i].facets.size(); ++j) {
                        if (are_equal_facets(facet, vertices, polyhedrons[i].facets[j], polyhedrons[i].vertices)) {

                            neigh.first  = i;
                            neigh.second = j;

                            return;
                        }
                    }
                }
            }

            bool are_equal_facets(const Polyhedron_facet &f1, const Polyhedron_vertices &v1, const Polyhedron_facet &f2, const Polyhedron_vertices &v2) const {

                if (f1.indices.size() != f2.indices.size()) return false;

                size_t count = 0;
                for (size_t i = 0; i < f1.indices.size(); ++i) {
                    for (size_t j = 0; j < f2.indices.size(); ++j) {

                        if (CGAL::abs(v1[f1.indices[i]].x() - v2[f2.indices[j]].x()) < m_tolerance && 
                            CGAL::abs(v1[f1.indices[i]].y() - v2[f2.indices[j]].y()) < m_tolerance && 
                            CGAL::abs(v1[f1.indices[i]].z() - v2[f2.indices[j]].z()) < m_tolerance) {
                            ++count; break;
                        }
                    }
                }
                return count == f1.indices.size();
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_WEIGHT_QUALITY_POLYHEDRON_ESTIMATOR_STEP_10_H