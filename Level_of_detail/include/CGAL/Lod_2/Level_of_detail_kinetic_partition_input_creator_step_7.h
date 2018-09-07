#ifndef CGAL_LEVEL_OF_DETAIL_KINETIC_PARTITION_INPUT_CREATOR_STEP_7_H
#define CGAL_LEVEL_OF_DETAIL_KINETIC_PARTITION_INPUT_CREATOR_STEP_7_H

// STL includes.
#include <vector>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>

// Jean Philippe includes.
#include <CGAL/Buildings/jean_philippe/defs_cgal.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_kinetic_partition_input_creator_step_7 {
            
        public:
            using Kernel    = InputKernel;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT        = typename Kernel::FT;
            using Point_2   = typename Kernel::Point_2;
            using Point_3   = typename Kernel::Point_3;
            using Vector_3  = typename Kernel::Vector_3;

            using Buildings_iterator = typename Buildings::iterator;

            using Floor_faces = typename Building::Floor_faces;
            
            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;
            
            using Roof  = typename Building::Roof;
            using Roofs = typename Building::Roofs;
            
            using JP_polygon  = typename Building::JP_polygon;
            using JP_polygons = typename Building::JP_polygons;

            using JP_FT      = JPTD::FT;
            using JP_point_3 = JPTD::CGAL_Point_3;

            using Points_2 = std::vector<Point_2>;
            using Points_3 = std::vector<Point_3>;

            using Polygon = Points_3;

            typename Kernel::Compute_squared_length_3   squared_length_3;
            typename Kernel::Compute_squared_distance_3 squared_distance_3;

            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Local_kernel = CGAL::Simple_cartesian<double>;
			using Local_ft     = Local_kernel::FT;
            using Point_2ft    = Local_kernel::Point_2;

            using Point_creator_2ft = Creator_uniform_2<Local_ft, Point_2ft>;
            using Points_2ft = std::vector<Point_2ft>;

            Level_of_detail_kinetic_partition_input_creator_step_7(const FT ground_height, Buildings &buildings) :
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_big_value(FT(100000000000000)),
            m_z_scale(FT(10)),
            m_up_scale(FT(3) / FT(2)),
            m_down_scale(FT(1) / FT(2)),
            m_disc_scale(FT(1) / FT(10)),
            m_num_points_in_disc(25),
            m_perturbation_axis_range(100),
            m_max_perturbation_angle_deg(3),
            m_scale_polygon(true),
            m_perturb_polygon_vertices(true), 
            m_perturb_polygon_support_plane(false),
            m_use_fixed_disc_radius(true),
            m_fixed_disc_radius(FT(1) / FT(1000)) {

                srand(time(NULL));
            }

            void scale(const bool new_state) {
                m_scale_polygon = new_state;
            }

            void perturb_vertices(const bool new_state) {
                m_perturb_polygon_vertices = new_state;
            }

            void perturb_support_plane(const bool new_state) {
                m_perturb_polygon_support_plane = new_state;
            }

            void create() const {
                
                if (m_buildings.size() == 0)
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;
                    
					if (building.is_valid && building.shapes.size() != 0) 
                        process_building(building);
                } 
            }

        private:
            const FT m_ground_height;
            Buildings &m_buildings;
            
            const FT m_big_value;
            const FT m_z_scale;
            const FT m_up_scale;
            const FT m_down_scale;
            const FT m_disc_scale;
            
            const size_t m_num_points_in_disc;
            const size_t m_perturbation_axis_range;
            const size_t m_max_perturbation_angle_deg;

            bool m_scale_polygon;
            bool m_perturb_polygon_vertices;
            bool m_perturb_polygon_support_plane;

            const bool m_use_fixed_disc_radius;
            const FT   m_fixed_disc_radius;

            void process_building(Building &building) const {

                JP_polygons &jp_polygons = building.jp_polygons;
                jp_polygons.clear();

                set_ground(building, jp_polygons);
                 set_walls(building, jp_polygons);
                 set_roofs(building, jp_polygons);
            }

            void set_ground(const Building &building, JP_polygons &jp_polygons) const {

                FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;

                const Floor_faces &faces = building.faces;
                for (size_t i = 0; i < faces.size(); ++i) {
                    
                    for (size_t j = 0; j < 3; ++j) {    
                        const Point_2 &p = faces[i]->vertex(j)->point();
                        
                        minx = CGAL::min(minx, p.x());
                        miny = CGAL::min(miny, p.y());

                        maxx = CGAL::max(maxx, p.x());
                        maxy = CGAL::max(maxy, p.y());
                    }
                }

                Polygon polygon(4);

                polygon[0] = Point_3(minx, miny, m_ground_height);
                polygon[1] = Point_3(maxx, miny, m_ground_height);
                polygon[2] = Point_3(maxx, maxy, m_ground_height);
                polygon[3] = Point_3(minx, maxy, m_ground_height);

                process_polygon(polygon, jp_polygons, m_up_scale, FT(1));
            }

            void process_polygon(Polygon &polygon, JP_polygons &jp_polygons, const FT scale, const FT z_extender) const {

                if (polygon.size() == 0) return;

                if (m_scale_polygon)                 scale_polygon(scale, z_extender, polygon);
                if (m_perturb_polygon_vertices)      perturb_polygon_vertices(polygon);
                if (m_perturb_polygon_support_plane) perturb_polygon_support_plane(polygon);
                
                JP_polygon jp_polygon(polygon.size());

                for (size_t i = 0; i < polygon.size(); ++i)
                    jp_polygon[i] = JP_point_3(JP_FT(polygon[i].x()), JP_FT(polygon[i].y()), JP_FT(polygon[i].z()));

                jp_polygons.push_back(jp_polygon);
            }

            void scale_polygon(const FT scale, const FT z_extender, Polygon &polygon) const {

                Point_3 b;
                compute_barycentre(polygon, b);

                for (size_t i = 0; i < polygon.size(); ++i) {
                    Point_3 &p = polygon[i];

                    const FT x = (p.x() - b.x()) * scale + b.x();
                    const FT y = (p.y() - b.y()) * scale + b.y();
                    const FT z = (p.z() - b.z()) * scale * z_extender + b.z();

                    p = Point_3(x, y, z);
                }
            }

            void compute_barycentre(const Polygon &polygon, Point_3 &b) const {

                CGAL_precondition(polygon.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < polygon.size(); ++i) {
                    const Point_3 &p = polygon[i];

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(polygon.size());
                y /= static_cast<FT>(polygon.size());
                z /= static_cast<FT>(polygon.size());

                b = Point_3(x, y, z);
            }

            void perturb_polygon_vertices(Polygon &polygon) const {

                Vector_3 source_normal;
                compute_source_normal(polygon, source_normal);

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                FT angle;
                Vector_3 axis;
                compute_angle_and_axis(source_normal, target_normal, angle, axis);

                Point_3 b;
                compute_barycentre(polygon, b);

                if (angle != FT(0)) rotate_polygon(b, angle, axis, polygon);
                
                perturb_horizontal_polygon(polygon);

                if (angle != FT(0)) rotate_polygon(b, -angle, axis, polygon);
            }

            void compute_source_normal(const Polygon &polygon, Vector_3 &normal) const {

                CGAL_precondition(polygon.size() >= 3);
                const Point_3 &p1 = polygon[0];
                const Point_3 &p2 = polygon[1];
                const Point_3 &p3 = polygon[2];

                const Vector_3 v1 = Vector_3(p1, p2);
                const Vector_3 v2 = Vector_3(p1, p3);

                normal = cross_product_3(v1, v2);
                normalize(normal);
            }

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void compute_target_normal(Vector_3 &normal) const {
                normal = Vector_3(FT(0), FT(0), FT(1));
            }

            void compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {
				
				const auto cross = cross_product_3(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const FT dot     = dot_product_3(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                if (angle == FT(0)) return;

                CGAL_precondition(length != FT(0));
                axis = cross / length;

                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle > half_pi) {
                    
                    angle = static_cast<FT>(CGAL_PI) - angle;
                    axis = -axis;
                }
			}

            void rotate_polygon(const Point_3 &b, const FT angle, const Vector_3 &axis, Polygon &polygon) const {

                Point_3 q;
                for (size_t i = 0; i < polygon.size(); ++i) {   
                    Point_3 &p = polygon[i];

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

            void perturb_horizontal_polygon(Polygon &polygon) const {

                const FT disc_radius = compute_disc_radius(polygon); 
                for (size_t i = 0; i < polygon.size(); ++i) {

                    Point_3 &p = polygon[i];
                    perturb_point_inside_disc(disc_radius, p);
                }
            }

            FT compute_disc_radius(const Polygon &polygon) const {

                if (m_use_fixed_disc_radius) 
                    return m_fixed_disc_radius;
                
                FT min_dist = m_big_value;

                for (size_t i = 0; i < polygon.size(); ++i) {
                    const size_t ip = (i + 1) % polygon.size();

                    const Point_3 &p1 = polygon[i];
                    const Point_3 &p2 = polygon[ip];

                    const FT dist = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_3(p1, p2))));
                    min_dist = CGAL::min(min_dist, dist);
                }

                return min_dist * m_disc_scale;
            }

            void perturb_point_inside_disc(const FT disc_radius, Point_3 &p) const {

                Points_2ft points_2ft;
                points_2ft.reserve(m_num_points_in_disc);

                Random_points_in_disc_2<Point_2ft, Point_creator_2ft> random_points_2ft(CGAL::to_double(disc_radius));
                CGAL::cpp11::copy_n(random_points_2ft, m_num_points_in_disc, std::back_inserter(points_2ft));

                const size_t rand_index = size_t_rand(m_num_points_in_disc - 1);
                const Point_2ft &q = points_2ft[rand_index];

                const FT qx = static_cast<FT>(q.x());
                const FT qy = static_cast<FT>(q.y());

                p = Point_3(p.x() + qx, p.y() + qy, p.z());
            }

            size_t size_t_rand(const size_t maxv) const {
                return static_cast<size_t>(rand() % maxv);
            }

            void perturb_polygon_support_plane(Polygon &polygon) const {

                const FT x = static_cast<FT>(size_t_rand(m_perturbation_axis_range));
                const FT y = static_cast<FT>(size_t_rand(m_perturbation_axis_range));
                const FT z = static_cast<FT>(size_t_rand(m_perturbation_axis_range));

                Vector_3 axis = Vector_3(x, y, z); 
                normalize(axis);

                const FT angle = static_cast<FT>(size_t_rand(m_max_perturbation_angle_deg)) * static_cast<FT>(CGAL_PI) / FT(180);
                if (angle != FT(0)) {

                    Point_3 b;
                    compute_barycentre(polygon, b);

                    rotate_polygon(b, angle, axis, polygon);
                }
            }

            void set_walls(Building &building, JP_polygons &jp_polygons) const {

                Walls &walls = building.walls;
                for (size_t i = 0; i < walls.size(); ++i) {
                    
                    Wall &wall = walls[i];
                    process_polygon(wall.boundary, jp_polygons, m_down_scale, m_z_scale);
                }
            }

            void set_roofs(Building &building, JP_polygons &jp_polygons) const {

                Roofs &roofs = building.roofs;
                for (size_t i = 0; i < roofs.size(); ++i) {
                    
                    if (!roofs[i].is_valid) 
                        continue;

                    Roof &roof = roofs[i];
                    process_polygon(roof.boundary, jp_polygons, m_down_scale, FT(1));
				}
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_KINETIC_PARTITION_INPUT_CREATOR_STEP_7_H