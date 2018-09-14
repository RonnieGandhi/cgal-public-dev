#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_15_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_15_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Lod_2/Facet_regions_merger_3/Level_of_detail_facet_regions_merger_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_coplanar_facets_merger_step_15 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using Buildings_iterator = typename Buildings::iterator;
            using Clean_facets       = typename Building::Clean_facets;
            using Output_regions     = typename Building::Output_regions;

            using Facet_regions_merger = Level_of_detail_facet_regions_merger_3<Kernel, Building>;
			
            Level_of_detail_coplanar_facets_merger_step_15(Buildings &buildings) :
            m_buildings(buildings)
            { }

            void merge() {
                
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

                Clean_facets &clean_facets = building.clean_facets;
                clean_facets.clear();

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(clean_facets);

                Output_regions &output_regions = building.output_regions;
                output_regions.clear();

                output_regions.push_back(clean_facets);
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_FACETS_MERGER_STEP_15_H