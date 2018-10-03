#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_14_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_14_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

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
            using Point_3  = typename Kernel::Point_3;
            using Vector_3 = typename Kernel::Vector_3;

            using Buildings_iterator = typename Buildings::iterator;
            
            using Clean_facet        = typename Building::Clean_facet;
            using Clean_facets       = typename Building::Clean_facets;

            using Output_regions = typename Building::Output_regions;

            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;

            using Roof  = typename Building::Roof;
            using Roofs = typename Building::Roofs;

            using Facet_regions_merger = Level_of_detail_facet_regions_merger_3<Kernel, Building>;

            using Log = CGAL::LOD::Mylog;

            typename Kernel::Compute_squared_length_3 squared_length_3;

            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;
			
            Level_of_detail_coplanar_facets_merger_step_14(Buildings &buildings, const FT ground_height) :
            m_buildings(buildings),
            m_use_original_facets(true),
            m_tolerance(FT(1) / FT(100000)),
            m_ground_height(ground_height)
            { }

            void merge() {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }

                create_walls_and_roofs();
                save_walls_and_roofs();
            }

            void use_merged_facets(const bool new_state) {
                m_use_original_facets = !new_state;
            }

        private:
            Buildings &m_buildings;
            bool m_use_original_facets;

            const FT m_tolerance;
            const FT m_ground_height;

            void process_building(Building &building) const {

                if (m_use_original_facets)
                    return;

                Clean_facets &clean_facets = building.clean_facets;
                clean_facets.clear();

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(clean_facets);

                Output_regions &output_regions = building.output_regions;
                output_regions.clear();

                output_regions.push_back(clean_facets);
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
                    if (!are_equal_points(facet_normal, zero)) {
                     
                        normalize(facet_normal);
                        return true;
                    }
                }
                return false;
			}

            template<class Point>
            bool are_equal_points(const Point &p, const Point &q) const {

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
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_14_H