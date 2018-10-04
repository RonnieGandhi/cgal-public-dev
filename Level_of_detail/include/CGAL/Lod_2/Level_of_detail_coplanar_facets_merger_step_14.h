#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_14_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_14_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Lod_2/Facet_regions_merger_3/Level_of_detail_facet_regions_merger_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_coplanar_facets_merger_step_14 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT       = typename Kernel::FT;
            using Point_2  = typename Kernel::Point_2;
            using Point_3  = typename Kernel::Point_3;
            using Vector_3 = typename Kernel::Vector_3;

            using Buildings_iterator = typename Buildings::iterator;
            
            using Clean_facet  = typename Building::Clean_facet;
            using Clean_facets = typename Building::Clean_facets;

            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;

            using Roof  = typename Building::Roof;
            using Roofs = typename Building::Roofs;

            using Boundary = typename Building::Boundary;

            using Facet_regions_merger = Level_of_detail_facet_regions_merger_3<Kernel, Building>;

            using Log = CGAL::LOD::Mylog;

            typename Kernel::Compute_squared_length_3 squared_length_3;

            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Index_pair       = std::pair<size_t, size_t>;
            using Point_with_index = std::pair<Point_3, Index_pair>;
            using Neighbours       = std::vector<Point_with_index>;

			using Point_map = CGAL::First_of_pair_property_map<Point_with_index>;

            using Search_traits_3 = CGAL::Search_traits_3<Kernel>;
            using Search_traits   = CGAL::Search_traits_adapter<Point_with_index, Point_map, Search_traits_3>;
			using Fuzzy_tree      = CGAL::Kd_tree<Search_traits>;
            using Fuzzy_sphere    = CGAL::Fuzzy_sphere<Search_traits>;
			
            Level_of_detail_coplanar_facets_merger_step_14(Buildings &buildings, const FT ground_height) :
            m_buildings(buildings),
            m_use_original_facets(true),
            m_tolerance(FT(1) / FT(100000)),
            m_ground_height(ground_height),
            m_scale(-FT(1)),
            m_proximity(-FT(1)),
            m_optimize_vertices(false)
            { }

            void merge() {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
                
                if (m_optimize_vertices) 
                    optimize_all_vertices();

                create_walls_and_roofs();
                save_walls_and_roofs();
            }

            void set_scale(const FT new_value) {
                m_scale = new_value;
            }

            void set_proximity(const FT new_value) {
                m_proximity = new_value;
            }

            void optimize_vertices(const bool new_value) {
                m_optimize_vertices = new_value;
            }

            void use_merged_facets(const bool new_state) {
                m_use_original_facets = !new_state;
            }

        private:
            Buildings &m_buildings;
            bool m_use_original_facets;

            const FT m_tolerance;
            const FT m_ground_height;
            
            FT m_scale;
            FT m_proximity;

            bool m_optimize_vertices;

            void process_building(Building &building) const {

                if (m_use_original_facets)
                    return;

                Clean_facets &clean_facets = building.clean_facets;
                clean_facets.clear();

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(clean_facets);
            }

            void create_walls_and_roofs() const {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building_walls_and_roofs(building);
                }
            }

            void process_building_walls_and_roofs(Building &building) const {

                create_walls(building);
                create_roofs(building);
            }

            void create_walls(Building &building) const {

				Walls &walls = building.walls;
                walls.clear();
                
                const Clean_facets &clean_facets = building.clean_facets;
				for (size_t i = 0; i < clean_facets.size(); ++i) {
				
					if (is_wall_facet(clean_facets[i]))
                        add_wall(clean_facets[i], walls);
				}
            }

            bool is_wall_facet(const Clean_facet &vertices) const {

				Vector_3 facet_normal;
				const bool success = set_facet_normal(vertices, facet_normal);

				if (!success) return false;

				Vector_3 ground_normal;
				set_ground_normal(ground_normal);

                const FT angle      = compute_angle(facet_normal, ground_normal);
                const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle));

                if (angle_diff < FT(5)) return true;
                return false;
            }

            bool set_facet_normal(const Clean_facet &vertices, Vector_3 &facet_normal) const {
                
                CGAL_precondition(vertices.size() >= 3);

                const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
                for (size_t i = 0; i < vertices.size(); ++i) {

                    const size_t ip  = (i + 1) % vertices.size();
                    const size_t ipp = (i + 2) % vertices.size();
                    
                    const Point_3 &p1 = vertices[i];
                    const Point_3 &p2 = vertices[ip];
                    const Point_3 &p3 = vertices[ipp];

                    const Vector_3 v1 = Vector_3(p2, p1);
                    const Vector_3 v2 = Vector_3(p2, p3);

                    facet_normal = cross_product_3(v1, v2);
                    if (!are_equal_points_3(facet_normal, zero)) {
                     
                        normalize(facet_normal);
                        return true;
                    }
                }
                return false;
			}

            template<class Point>
            bool are_equal_points_3(const Point &p, const Point &q) const {

                const FT eps = m_tolerance;
                return (CGAL::abs(p.x() - q.x()) < eps) && (CGAL::abs(p.y() - q.y()) < eps) && (CGAL::abs(p.z() - q.z()) < eps);
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

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void add_wall(const Clean_facet &clean_facet, Walls &walls) const {

                Wall wall;
                wall.boundary = clean_facet;
                walls.push_back(wall);
            }

            void create_roofs(Building &building) const {

                Roofs &roofs = building.roofs;
                roofs.clear();
                
                const Clean_facets &clean_facets = building.clean_facets;
				for (size_t i = 0; i < clean_facets.size(); ++i) {
				
					if (!is_wall_facet(clean_facets[i]) && !is_ground_facet(clean_facets[i]))
                        add_roof(clean_facets[i], roofs);
				}
            }

			bool is_ground_facet(const Clean_facet &clean_facet) const {
				
				FT average_height = FT(0);
				for (size_t i = 0; i < clean_facet.size(); ++i) {
					
					const Point_3 &p = clean_facet[i];
					average_height += p.z();
				}
				
				average_height /= static_cast<FT>(clean_facet.size());
				return CGAL::abs(average_height - m_ground_height) < m_tolerance * FT(100);
			}

            void add_roof(const Clean_facet &clean_facet, Roofs &roofs) const {

                Roof roof;
                roof.boundary = clean_facet;
                roofs.push_back(roof);
            }

            void save_walls_and_roofs() const {

                Log exporter; 
                exporter.save_building_walls(m_buildings, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "15_walls", true);
                exporter.save_building_roofs_without_faces(m_buildings, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "16_roofs", true);
            }

            void optimize_all_vertices() {

                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        optimize_building_vertices(building);
                }
            }

            void optimize_building_vertices(Building &building) const {

                const Clean_facets &clean_facets = building.clean_facets;
                
                Fuzzy_tree tree;
                create_tree(clean_facets, tree);

				for (size_t i = 0; i < clean_facets.size(); ++i)
                    process_clean_facet(tree, clean_facets[i], building);
            }

            void create_tree(const Clean_facets &clean_facets, Fuzzy_tree &tree) const {

                tree.clear();
                for (size_t i = 0; i < clean_facets.size(); ++i)
                    for (size_t j = 0; j < clean_facets[i].size(); ++j)
                        tree.insert(std::make_pair(clean_facets[i][j], std::make_pair(i, j)));
            }

            void process_clean_facet(const Fuzzy_tree &tree, const Clean_facet &clean_facet, Building &building) const {

                for (size_t i = 0; i < clean_facet.size(); ++i)
                    process_clean_facet_vertex(tree, clean_facet[i], building);
            }

            void process_clean_facet_vertex(const Fuzzy_tree &tree, const Point_3 &q, Building &building) const {

                const Fuzzy_sphere sphere = Fuzzy_sphere(q, m_scale);
                
                Neighbours neighbours;
                tree.search(std::back_inserter(neighbours), sphere);

                if (is_corner_vertex(q, building)) process_corner_vertex(q, neighbours, building.clean_facets);
                else process_standard_vertex(q, neighbours, building);
            }

            bool is_corner_vertex(const Point_3 &q, const Building &building) const {

                const Boundary &boundary = building.boundaries[0];
                for (size_t i = 0; i < boundary.size(); ++i) {
                    
                    const Point_2 &p = boundary[i]->point();
                    if (are_equal_points_2(p, Point_2(q.x(), q.y()))) return true;
                }
                return false;
            }

            template<class Point>
            bool are_equal_points_2(const Point &p, const Point &q) const {

                const FT eps = m_proximity;
                return (CGAL::abs(p.x() - q.x()) < eps) && (CGAL::abs(p.y() - q.y()) < eps);
            }

            void process_corner_vertex(const Point_3 &q, const Neighbours &neighbours, Clean_facets &clean_facets) const {

                for (size_t i = 0; i < neighbours.size(); ++i) {
                    
                    const auto &pair = neighbours[i].second;
                    clean_facets[pair.first][pair.second] = q;
                }
            }

            void process_standard_vertex(const Point_3 &q, const Neighbours &neighbours, Building &building) const {

                for (size_t i = 0; i < neighbours.size(); ++i) {
                    const auto &pair = neighbours[i].second;

                    if (is_corner_vertex(building.clean_facets[pair.first][pair.second], building))
                        return;
                }

                FT x = q.x(), y = q.y(), z = q.z();

                for (size_t i = 0; i < neighbours.size(); ++i) {
                    const auto &pair = neighbours[i].second;

                    x += building.clean_facets[pair.first][pair.second].x();
                    y += building.clean_facets[pair.first][pair.second].y();
                    z += building.clean_facets[pair.first][pair.second].z();
                }

                x /= static_cast<FT>(neighbours.size() + 1);
                y /= static_cast<FT>(neighbours.size() + 1);
                z /= static_cast<FT>(neighbours.size() + 1);

                const Point_3 res = Point_3(x, y, z);

                for (size_t i = 0; i < neighbours.size(); ++i) {
                    
                    const auto &pair = neighbours[i].second;
                    building.clean_facets[pair.first][pair.second] = res;
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_14_H