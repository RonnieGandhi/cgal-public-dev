#ifndef CGAL_LEVEL_OF_DETAIL_POLYHEDRON_IN_OUT_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_POLYHEDRON_IN_OUT_ESTIMATOR_H 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding>
		class Level_of_detail_polyhedron_in_out_estimator {

		public:
            using Kernel   = InputKernel;
			using Input    = InputData;
			using Building = InputBuilding;

			using FT 	     = typename Kernel::FT;
			using Line_2     = typename Kernel::Line_2;
            using Line_3     = typename Kernel::Line_3;
			using Point_2    = typename Kernel::Point_2;
			using Point_3    = typename Kernel::Point_3;
            using Plane_3    = typename Kernel::Plane_3;
            using Vector_3   = typename Kernel::Vector_3;
			using Triangle_2 = typename Kernel::Triangle_2;

            using Index       = typename Building::Index;
            using Indices     = typename Building::Indices;
			using Polyhedron  = typename Building::Polyhedron;
			using Polyhedrons = typename Building::Polyhedrons;
			using Floor_faces = typename Building::Floor_faces;

			using Vertices = typename Polyhedron::Vertices;
            using Facet    = typename Polyhedron::Facet;
            using Facets   = typename Polyhedron::Facets;

			using Pair       = CGAL::cpp11::array<FT, 2>;
            using Statistics = std::pair<FT, FT>;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;
            typename Kernel::Compute_squared_distance_3 squared_distance_3;
            
            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
            typename Kernel::Construct_cross_product_vector_3 cross_product_3;
            typename Kernel::Compute_squared_length_3         squared_length_3;

            using Points_3    = std::vector<Point_3>;
            using Polygon     = CGAL::Polygon_2<Kernel>;
            using Coordinates = std::vector<FT>;

            using Mean_value = CGAL::Barycentric_coordinates::Mean_value_2<Kernel>;
            using Mean_value_coordinates = CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel>;

            using Intersect = typename Kernel::Intersect_3;

			Level_of_detail_polyhedron_in_out_estimator(const Input &input, const FT ground_height) : 
			m_input(input),
			m_ground_height(ground_height),
			m_big_value(FT(100000000000000)),
			m_distance_tolerance(FT(3) / FT(10)),
			m_bc_tolerance_top(FT(6) / FT(5)),
            m_bc_tolerance_bottom(-FT(1) / FT(5)),
            m_angle_threshold(FT(10)),
            m_edge_length_tolerance(FT(1) / FT(1000)),
            m_height_offset(FT(1) / FT(8)) { }

			void estimate(Building &building) const {
				compute_building_maximum_height(building);

				Polyhedrons &polyhedrons = building.polyhedrons;
                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    
					Polyhedron &polyhedron = polyhedrons[i];
					polyhedron.out = estimate_out_value(building, polyhedron);

					CGAL_precondition(polyhedron.out >= FT(0) && polyhedron.out <= FT(1));
					polyhedron.in = FT(1) - polyhedron.out;
                }
			}

		private:
			const Input &m_input;

			const FT m_ground_height;
			const FT m_big_value;
			const FT m_distance_tolerance;
			const FT m_bc_tolerance_top;
            const FT m_bc_tolerance_bottom;
            const FT m_angle_threshold;
            const FT m_edge_length_tolerance;
            const FT m_height_offset;

            void compute_building_maximum_height(Building &building) const {

                const Indices &interior_indices = building.interior_indices;
                CGAL_precondition(interior_indices.size() > 0);

                FT max_height = -m_big_value;
				for (size_t i = 0; i < interior_indices.size(); ++i) {
                    
                    const Index point_index = interior_indices[i];
                    const Point_3 &p = m_input.point(point_index);

                    max_height = CGAL::max(max_height, p.z());
                }

                building.max_height = max_height;
                CGAL_postcondition(max_height >= FT(0));
            }

			FT estimate_out_value(const Building &building, const Polyhedron &polyhedron) const {

				Point_3 barycentre;
                compute_polyhedron_barycentre(polyhedron, barycentre);

                if (is_above_building_max_height(barycentre, building.max_height)) return FT(1);
                if (is_below_ground(barycentre, m_ground_height)) 				   return FT(1);
                if (is_out_of_building(barycentre, building)) 					   return FT(4) / FT(5);
                if (has_vertices_outside(polyhedron, building))                    return FT(3) / FT(5);

                const Statistics stats = get_statistics(polyhedron, building);
                
                const FT in  = stats.first;
                const FT out = stats.second;

                if (in == FT(0) && out == FT(0)) return 0;
                return out / (in + out);
			}

			void compute_polyhedron_barycentre(const Polyhedron &polyhedron, Point_3 &barycentre) const {

                const Vertices &vertices  = polyhedron.vertices;
                const size_t num_vertices = vertices.size();

                CGAL_precondition(num_vertices > 0);

                FT x = FT(0), y = FT(0), z = FT(0);
                for (size_t i = 0; i < num_vertices; ++i) {

                    x += vertices[i].x();
                    y += vertices[i].y();
                    z += vertices[i].z();
                }

                x /= static_cast<FT>(num_vertices);
                y /= static_cast<FT>(num_vertices);
                z /= static_cast<FT>(num_vertices);

                barycentre = Point_3(x, y, z);
            }

			bool is_above_building_max_height(const Point_3 &query, const FT building_max_height) const {
                return query.z() > building_max_height;
            }

			bool is_below_ground(const Point_3 &query, const FT ground_height) const {
                return query.z() < ground_height;
            }

			bool is_out_of_building(const Point_3 &query, const Building &building) const {
                
                const Point_2 p = Point_2(query.x(), query.y());

                const Floor_faces &floor_faces = building.faces;
                for (size_t i = 0; i < floor_faces.size(); ++i) {
                        
                    const Point_2 &p1 = floor_faces[i]->vertex(0)->point();
                    const Point_2 &p2 = floor_faces[i]->vertex(1)->point();
                    const Point_2 &p3 = floor_faces[i]->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (is_within_triangle(p, triangle)) return false;
                }
                return true;
            }

			bool is_within_triangle(const Point_2 &query, const Triangle_2 &triangle) const {
                
                if (triangle.has_on_bounded_side(query) || triangle.has_on_boundary(query)) 
                    return true;
                
                for (size_t i = 0; i < 3; ++i) {
                    const size_t ip = (i + 1) % 3;

                    const Point_2 &p1 = triangle.vertex(i);
                    const Point_2 &p2 = triangle.vertex(ip);

                    const Line_2 line = Line_2(p1, p2);

                    const Point_2 projected   = line.projection(query);
                    const FT squared_distance = squared_distance_2(query, projected);

                    const Pair pair = BC::compute_segment_coordinates_2(p1, p2, projected, Kernel());

                    const FT squared_tolerance = m_distance_tolerance * m_distance_tolerance;
                    
                    const FT epst = m_bc_tolerance_top;
                    const FT epsb = m_bc_tolerance_bottom;

                    if (pair[0] > epsb && pair[1] > epsb && pair[0] < epst && pair[1] < epst && squared_distance < squared_tolerance) return true;
                }
                return false;
            }

            bool has_vertices_outside(const Polyhedron &polyhedron, const Building &building) const {

                size_t count = 0;
                const Vertices &vertices = polyhedron.vertices;
                for (size_t i = 0; i < vertices.size(); ++i) {

                    const Point_3 &vertex = vertices[i];
                    const bool is_out = is_out_of_building(vertex, building);

                    if (is_out) ++count;
                    if (is_out && count > 0) return true;
                }
                return false;
            }

            Statistics get_statistics(const Polyhedron &polyhedron, const Building &building) const {

                size_t in = 0, out = 0;
                Indices tmp_indices;

                process_facets(polyhedron, building, tmp_indices, in, out);
                process_middle_plane(polyhedron, tmp_indices, in, out);

                return std::make_pair(static_cast<FT>(in), static_cast<FT>(out));
            }

            void process_facets(const Polyhedron &polyhedron, const Building &building, Indices &tmp_indices, size_t &in, size_t &out) const {

                const Vertices &vertices = polyhedron.vertices;
                const Facets   &facets   = polyhedron.facets;

                for (size_t i = 0; i < facets.size(); ++i) {    
                    const Facet &facet = facets[i];

                    if (is_vertical_facet(vertices, facet)) continue;
                    process_facet(vertices, facet, building, tmp_indices, in, out);
                }
            }

            bool is_vertical_facet(const Vertices &vertices, const Facet &facet) const {

				Vector_3 facet_normal;
				const bool success = set_facet_normal(vertices, facet, facet_normal);

                if (!success) return true;

				Vector_3 ground_normal;
				set_ground_normal(ground_normal);

                const FT angle      = compute_angle(facet_normal, ground_normal);
                const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle));

                if (angle_diff < m_angle_threshold) return true;
                return false;
            }

            bool set_facet_normal(const Vertices &vertices, const Facet &facet, Vector_3 &facet_normal) const {

                CGAL_precondition(facet.indices.size() >= 3);

                Points_3 points;
                get_points(vertices, facet, points);

                if (points.size() < 3) return false;

                const Point_3 &p1 = points[0];
                const Point_3 &p2 = points[1];
                const Point_3 &p3 = points[2];

                const Vector_3 v1 = Vector_3(p1, p2);
                const Vector_3 v2 = Vector_3(p1, p3);

                facet_normal = cross_product_3(v1, v2);
                normalize(facet_normal);

                return true;
			}

            void get_points(const Vertices &vertices, const Facet &facet, Points_3 &points) const {

                points.clear();
                const size_t num_indices = facet.indices.size();

                for (size_t i = 0; i < num_indices; ++i) {
                    const size_t ip = (i + 1) % num_indices;

                    const Point_3 &p1 = vertices[facet.indices[i]];
                    const Point_3 &p2 = vertices[facet.indices[ip]];

                    const FT edge_length = static_cast<FT>(std::sqrt(CGAL::to_double(squared_distance_3(p1, p2))));

                    if (edge_length < m_edge_length_tolerance) continue;
                    else points.push_back(p1);
                }
            }

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void set_ground_normal(Vector_3 &ground_normal) const {
				ground_normal = Vector_3(FT(0), FT(0), FT(1));
			}

            FT compute_angle(const Vector_3 &m, const Vector_3 &n) const {
				
				const Vector_3 cross  = cross_product_3(m, n);
				const FT       length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const FT       dot    = dot_product_3(m, n);

				FT angle_rad = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                
                const FT      pi = static_cast<FT>(CGAL_PI);
                const FT half_pi = pi / FT(2);

                if (angle_rad > half_pi) angle_rad = pi - angle_rad;
                const FT angle_deg = angle_rad * FT(180) / pi;
                
                return angle_deg;
			}

            void process_facet(const Vertices &vertices, const Facet &facet, const Building &building, Indices &tmp_indices, size_t &in, size_t &out) const {

                Polygon polygon;
                const bool success = create_polygon(vertices, facet, polygon);

                if (!success) return;
                Mean_value_coordinates mean_value_coordinates(polygon.vertices_begin(), polygon.vertices_end());

                const auto &shapes = building.shapes;
                for (size_t i = 0; i < shapes.size(); ++i) {

                    const Indices &indices = shapes[i];
                    for (size_t j = 0; j < indices.size(); ++j) {
                    
                        const Index point_index = indices[j];
                        const Point_3 &p = m_input.point(point_index);

                        const FT height = intersect_facet(p, vertices, facet, mean_value_coordinates);
                        
                        if (height == m_big_value) continue;
                        tmp_indices.push_back(point_index);

                        if (is_inside_building(height, p.z())) ++in;
                        else ++out;
                    }
                }
            }

            bool create_polygon(const Vertices &vertices, const Facet &facet, Polygon &polygon) const {

                Points_3 points;
                get_points(vertices, facet, points);

                if (points.size() < 3) return false;

                polygon.clear();
                for (size_t i = 0; i < points.size(); ++i) {
                    
                    const Point_3 &p = points[i];
                    polygon.push_back(Point_2(p.x(), p.y()));
                }
                
                if (polygon.is_clockwise_oriented()) 
                    polygon.reverse_orientation();

                return polygon.is_simple();
            }

            FT intersect_facet(const Point_3 &p, const Vertices &vertices, const Facet &facet, Mean_value_coordinates &mean_value_coordinates) const {

                if (is_inside_polygon(p, mean_value_coordinates)) 
                    return intersect_line_and_plane(vertices, facet, p);
                
                return m_big_value;
            }

            bool is_inside_polygon(const Point_3 &p, Mean_value_coordinates &mean_value_coordinates) const {

                Coordinates coordinates;
                const Point_2 query = Point_2(p.x(), p.y());
                
                mean_value_coordinates(query, std::back_inserter(coordinates));
                CGAL_precondition(coordinates.size() >= 3);

                for (size_t i = 0 ; i < coordinates.size(); ++i)
                    if (coordinates[i] <= FT(0) || coordinates[i] >= FT(1)) return false;
                return true;
            }

            FT intersect_line_and_plane(const Vertices &vertices, const Facet &facet, const Point_3 &p) const {

                Line_3 line;
                bool success = create_line(p, line);

                if (!success) return m_big_value;

                Plane_3 plane;
                success = create_plane(vertices, facet, plane);

                if (!success) return m_big_value;

                return intersect(line, plane);
            }

            bool create_line(const Point_3 &p, Line_3 &line) const {
                
                const Point_3 p1 = Point_3(p.x(), p.y(), m_ground_height - FT(10));
                const Point_3 p2 = Point_3(p.x(), p.y(), m_ground_height + FT(10));

                line = Line_3(p1, p2);
                return true;
            }

            bool create_plane(const Vertices &vertices, const Facet &facet, Plane_3 &plane) const {

                CGAL_precondition(facet.indices.size() >= 3);

                Points_3 points;
                get_points(vertices, facet, points);

                if (points.size() < 3) return false;

                const Point_3 &p1 = points[0];
                const Point_3 &p2 = points[1];
                const Point_3 &p3 = points[2];

                plane = Plane_3(p1, p2, p3);
                return true;
            }

            FT intersect(const Line_3 &line, const Plane_3 &plane) const {

				typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result = intersection(line, plane);
                
                if (result) {

                    if (const Line_3 *tmp = boost::get<Line_3>(&*result)) return m_big_value;
                    else {
                        const Point_3 *point = boost::get<Point_3>(&*result);
				        return (*point).z();
                    }
                }
                return m_big_value;
            }

            bool is_inside_building(const FT current_height, const FT real_height) const {
                return current_height > m_ground_height - m_height_offset && current_height < real_height + m_height_offset;
            }

            void process_middle_plane(const Polyhedron &polyhedron, const Indices &indices, size_t &in, size_t &out) const {

                Point_3 barycentre;
                compute_polyhedron_barycentre(polyhedron, barycentre);

                Plane_3 middle_plane(barycentre, Vector_3(FT(0), FT(0), FT(1)));

                Line_3 line;
                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    const Index point_index = indices[i];
                    const Point_3 &p = m_input.point(point_index);

                    const bool success = create_line(p, line);
                    if (!success) continue;

                    const FT height = intersect(line, middle_plane);

                    if (is_inside_building(height, p.z())) ++in;
                    else ++out;
                }
            }
		};

	} // CGAL

} // LOD

#endif // CGAL_LEVEL_OF_DETAIL_POLYHEDRON_IN_OUT_ESTIMATOR_H