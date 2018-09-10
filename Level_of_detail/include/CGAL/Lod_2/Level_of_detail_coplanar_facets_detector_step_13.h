#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_DETECTOR_STEP_13_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_DETECTOR_STEP_13_H

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
		class Level_of_detail_coplanar_facets_detector_step_13 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using Buildings_iterator          = typename Buildings::iterator;
            using Facets_based_region_growing = CGAL::LOD::Level_of_detail_facets_based_region_growing_3<Kernel>;
            
            using Input_facets   = typename Building::Clean_facets;
            using Output_regions = typename Facets_based_region_growing::Output_regions;
			
            Level_of_detail_coplanar_facets_detector_step_13(Buildings &buildings) :
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

                const Input_facets &input_facets = building.clean_facets;
                Output_regions &output_regions   = building.output_regions;
                
                if (input_facets.size() == 0) return;

                Facets_based_region_growing region_growing(input_facets);
                region_growing.create_regions(output_regions);
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_DETECTOR_STEP_13_H