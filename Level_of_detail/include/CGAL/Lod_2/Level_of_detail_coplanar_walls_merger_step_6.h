#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_MERGER_STEP_6_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_MERGER_STEP_6_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Lod_2/Facet_regions_merger_3/Level_of_detail_facet_regions_merger_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_coplanar_walls_merger_step_6 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using Buildings_iterator = typename Buildings::iterator;

            using Walls         = typename Building::Walls;
            using Output_facet  = typename Building::Region_facet;
            using Output_facets = std::vector<Output_facet>;

            using Facet_regions_merger = Level_of_detail_facet_regions_merger_3<Kernel, InputBuilding>;
			
            Level_of_detail_coplanar_walls_merger_step_6(Buildings &buildings) :
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

                Walls &walls = building.walls;
                Output_facets output_facets;

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(output_facets);

                walls.clear();
                walls.resize(output_facets.size());

                for (size_t i = 0; i < output_facets.size(); ++i)
                    walls[i].boundary = output_facets[i];
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_MERGER_STEP_6_H