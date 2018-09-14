#ifndef CGAL_LEVEL_OF_DETAIL_ROOFS_REGULARIZER_STEP_4_H
#define CGAL_LEVEL_OF_DETAIL_ROOFS_REGULARIZER_STEP_4_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_roofs_regularizer_step_4 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using Buildings_iterator = typename Buildings::iterator;

            Level_of_detail_roofs_regularizer_step_4(Buildings &buildings) :
            m_buildings(buildings)
            { }

            void regularize() {
                
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

            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ROOFS_REGULARIZER_STEP_4_H