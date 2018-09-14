#ifndef CGAL_LEVEL_OF_DETAIL_WALLS_ESTIMATOR_STEP_5_H
#define CGAL_LEVEL_OF_DETAIL_WALLS_ESTIMATOR_STEP_5_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_walls_estimator_step_5 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using Buildings_iterator = typename Buildings::iterator;

            using Walls = typename Building::Walls;
            using Wall  = typename Building::Wall;

            using Wall_boundary     = typename Wall::Wall_boundary;
            using Building_boundary = typename Building::Boundary;

            Level_of_detail_walls_estimator_step_5(const FT ground_height, Buildings &buildings) :
            m_ground_height(ground_height),
            m_buildings(buildings)
            { }

            void estimate() {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
            }

        private:
            Buildings &m_buildings;
            const FT m_ground_height;

            void process_building(Building &building) const {
                
                Walls &walls = building.walls;
                walls.clear();

                create_walls(building, walls);
            }

            void create_walls(const Building &building, Walls &walls) const {

                const Building_boundary &boundary = building.boundaries[0];
                for (size_t i = 0; i < boundary.size(); i += 2) {
					
                    const size_t ip = i + 1;
					CGAL_precondition(ip < boundary.size());

                    const Point_2 &source =  boundary[i]->point();
                    const Point_2 &target = boundary[ip]->point();

                    create_wall(source, target, building.height, walls);
                }
            }

            void create_wall(const Point_2 &source, const Point_2 &target, const FT building_height, Walls &walls) const {

                Wall new_wall;
                Wall_boundary &boundary = new_wall.boundary;
                
                boundary.clear();
                boundary.resize(4);

                boundary[0] = Point_3(source.x(), source.y(), m_ground_height);
                boundary[1] = Point_3(target.x(), target.y(), m_ground_height);
                boundary[2] = Point_3(target.x(), target.y(), m_ground_height + building_height);
                boundary[3] = Point_3(source.x(), source.y(), m_ground_height + building_height);

                walls.push_back(new_wall);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_WALLS_ESTIMATOR_STEP_5_H