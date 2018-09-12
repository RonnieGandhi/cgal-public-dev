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
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

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

            using FT       = typename Kernel::FT;
            using Point_2  = typename Kernel::Point_2;
            using Point_3  = typename Kernel::Point_3;
            using Plane_3  = typename Kernel::Plane_3;
            using Vector_3 = typename Kernel::Vector_3;
            
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

            using Points_2 = std::vector<Point_2>;
            using Points_3 = std::vector<Point_3>;

            using Polygon = CGAL::Polygon_2<Kernel>;

            typename Kernel::Compute_squared_length_3   squared_length_3;
            typename Kernel::Compute_squared_distance_3 squared_distance_3;
            
            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Coordinates = std::vector<FT>;

            using Mean_value = CGAL::Barycentric_coordinates::Mean_value_2<Kernel>;
            using Mean_value_coordinates = CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel>;

            using Log = CGAL::LOD::Mylog;

            Level_of_detail_weight_quality_polyhedron_estimator_step_10(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height),
            m_tolerance(FT(1) / FT(100000)),
            m_edge_length_tolerance(FT(1) / FT(1000)),
            m_distance_threshold(FT(3))
            { }

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

            void process_building(Building &building) {

                Graphcut_facets &gc_facets = building.graphcut_facets;
                create_graphcut_facets(building, gc_facets);
            }

            void create_graphcut_facets(const Building &building, Graphcut_facets &gc_facets) const {
                
                gc_facets.clear();
                
                const Indices &interior_indices = building.interior_indices;
                const Polyhedrons &polyhedrons  = building.polyhedrons;

                for (int i = 0; i < polyhedrons.size(); ++i)
                    add_graphcut_facets(interior_indices, polyhedrons, i, gc_facets);
            }

            void add_graphcut_facets(const Indices &interior_indices, const Polyhedrons &polyhedrons, const int poly_index, Graphcut_facets &gc_facets) const {

                for (int i = 0; i < polyhedrons[poly_index].facets.size(); ++i)
                    process_polyhedron_facet(interior_indices, polyhedrons, poly_index, i, gc_facets);
            }

            void process_polyhedron_facet(const Indices &interior_indices, const Polyhedrons &polyhedrons, const int poly_index, const int facet_index, Graphcut_facets &gc_facets) const {

                if (was_already_added(poly_index, facet_index, gc_facets))
                    return;

                std::cout << poly_index << " " << facet_index << std::endl;

                if (poly_index == 2 && facet_index == 1) {

                    const Polyhedron &polyhedron = polyhedrons[poly_index];
                    
                    const Polyhedron_facets   &facets   = polyhedron.facets;
                    const Polyhedron_vertices &vertices = polyhedron.vertices;

                    Graphcut_facet gc_facet;
                    gc_facet.data = std::make_pair(poly_index, facet_index);

                    find_neighbours(polyhedrons, facets, vertices, poly_index, facet_index, gc_facet.neighbours);
                    
                    compute_weight(facets[facet_index], vertices, gc_facet);
                    compute_quality(interior_indices, facets[facet_index], vertices, gc_facet);

                    if (gc_facet.quality == FT(0)) gc_facets.push_back(gc_facet);
                    exit(0);
                }
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

            void compute_weight(const Polyhedron_facet &, const Polyhedron_vertices &, Graphcut_facet &gc_facet) const {
                gc_facet.weight = FT(1);
            }

            void compute_quality(const Indices &interior_indices, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Graphcut_facet &gc_facet) const {

                const Graphcut_facet_data &neigh_1 = gc_facet.neighbours.first;
                const Graphcut_facet_data &neigh_2 = gc_facet.neighbours.second;

                if (neigh_1.first < 0 || neigh_1.second < 0) {
                    
                    gc_facet.quality = FT(1);
                    return;
                }

                if (neigh_2.first < 0 || neigh_2.second < 0) {
                
                    gc_facet.quality = FT(1);
                    return;    
                }

                Indices closest;
                find_closest_points(interior_indices, facet, vertices, closest);

                if (closest.size() == 0) {
                    
                    gc_facet.quality = FT(1);
                    return;
                }

                Indices internal;
                find_internal_points(closest, facet, vertices, internal);

                if (internal.size() < 100) {
                    
                    gc_facet.quality = FT(1);
                    return;
                }

                // std::cout << internal.size() << std::endl;

                gc_facet.quality = FT(0);
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

            void find_internal_points(const Indices &closest, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Indices &internal) const {

                internal.clear();
                if (closest.size() == 0) return;

                Vector_3 source_normal;
                bool success = compute_source_normal(facet, vertices, source_normal);

                if (!success) return;

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                if (source_normal == -target_normal) source_normal = target_normal;

                FT angle; Vector_3 axis;
                success = compute_angle_and_axis(source_normal, target_normal, angle, axis);

                Point_3 b1;
                compute_barycentre(facet, vertices, b1);

                Point_3 b2;
                compute_barycentre(closest, b2);

                const Point_3 b = Point_3((b1.x() + b2.x()) / FT(2), (b1.y() + b2.y()) / FT(2), (b1.z() + b2.z()) / FT(2));

                if (!success) return;
                    
                Polygon polygon;
                success = create_polygon(facet, vertices, b1, angle, axis, polygon);

                if (!success) return;

                Points_2 queries;
                create_queries(closest, b2, angle, axis, queries);

                Log logger;
                logger.export_polygon("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_polygon", facet, vertices);
                logger.export_polygon("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_polygon_rotated", polygon);
                
                Points_3 closest_points(closest.size());
                for (size_t i = 0; i < closest.size(); ++i)
                    closest_points[i] = m_input.point(closest[i]);

                logger.export_3d_points("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_points", closest_points);
                logger.export_points("tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "debug_points_rotated", queries);

                find_points_inside_polygon(polygon, queries, closest, internal);
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

				const Vector_3 cross  = cross_product_3(m, n);
				const FT       length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const FT       dot    = dot_product_3(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                
                if (angle  == FT(0)) return true;
                if (length == FT(0)) return false;

                axis = cross / length;

                /*
                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle > half_pi) {
                    
                    angle = static_cast<FT>(CGAL_PI) - angle;
                    axis = -axis;
                } */

                return true;
			}

            void compute_barycentre(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Point_3 &b) const {
                
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

            void compute_barycentre(const Indices &closest, Point_3 &b) const {
                
                CGAL_precondition(closest.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < closest.size(); ++i) {
                    const Point_3 &p = m_input.point(closest[i]);

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(closest.size());
                y /= static_cast<FT>(closest.size());
                z /= static_cast<FT>(closest.size());

                b = Point_3(x, y, z);
            }

            bool create_polygon(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, const Point_3 &b, const FT angle, const Vector_3 &axis, Polygon &polygon) const {

                polygon.clear();

                for (size_t i = 0; i < facet.indices.size(); ++i) {
                    const Point_3 &p = vertices[facet.indices[i]];

                    if (angle != 0) {

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

                    if (angle != 0) {

                        Point_3 q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                        rotate_point(angle, axis, q);
                        queries.push_back(Point_2(q.x() + b.x(), q.y() + b.y()));

                    } else queries.push_back(Point_2(p.x(), p.y()));
                }
            }

            void find_points_inside_polygon(const Polygon &polygon, const Points_2 &queries, const Indices &closest, Indices &internal) const {
                
                CGAL_precondition(queries.size() == closest.size());

                internal.clear();
                Mean_value_coordinates mean_value_coordinates(polygon.vertices_begin(), polygon.vertices_end());

                for (size_t i = 0; i < queries.size(); ++i) {
                    const Point_2 &query = queries[i];

                    Coordinates coordinates;
                    mean_value_coordinates(query, std::back_inserter(coordinates));

                    if (is_inside_polygon(coordinates))
                        internal.push_back(closest[i]);
                }
            }

            bool is_inside_polygon(const Coordinates &coordinates) const {
                
                for (size_t i = 0 ; i < coordinates.size(); ++i)
                    if (coordinates[i] <= FT(0) || coordinates[i] >= FT(1)) return false;
                return true;
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