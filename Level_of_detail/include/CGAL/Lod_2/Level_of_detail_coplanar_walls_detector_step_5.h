#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_DETECTOR_STEP_5_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_DETECTOR_STEP_5_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Lod_2/Facets_based_region_growing_3/Level_of_detail_facets_based_region_growing_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_coplanar_walls_detector_step_5 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using Buildings_iterator = typename Buildings::iterator;

            using Walls = typename Building::Walls;

            using Facets_based_region_growing = CGAL::LOD::Level_of_detail_facets_based_region_growing_3<Kernel>;
            
            using Input_facets   = typename Facets_based_region_growing::Input_facets;
            using Output_regions = typename Facets_based_region_growing::Output_regions;
			
            Level_of_detail_coplanar_walls_detector_step_5(Buildings &buildings) :
            m_buildings(buildings)
            { }

            void detect() {
                
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

            void process_building(Building &building) const {

                Walls &walls                   = building.walls;
                Output_regions &output_regions = building.output_regions;
                
                Input_facets input_facets;
                create_input_facets(walls, input_facets);

                Facets_based_region_growing region_growing(input_facets);
                region_growing.create_regions(output_regions);
            }

            void create_input_facets(const Walls &walls, Input_facets &input_facets) const {
                
                input_facets.clear();
                input_facets.resize(walls.size());

                for (size_t i = 0; i < walls.size(); ++i)
                    input_facets[i] = walls[i].boundary;
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_DETECTOR_STEP_5_H