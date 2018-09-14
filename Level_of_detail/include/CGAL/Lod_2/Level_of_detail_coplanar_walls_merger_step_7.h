#ifndef CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_MERGER_STEP_7_H
#define CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_MERGER_STEP_7_H

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
		class Level_of_detail_coplanar_walls_merger_step_7 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT         = typename Kernel::FT;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_3 = typename Kernel::Triangle_3;

            using Buildings_iterator = typename Buildings::iterator;

            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;

            using Boundary = typename Wall::Wall_boundary;

            using Output_facet  = typename Building::Region_facet;
            using Output_facets = std::vector<Output_facet>;

            using Facet_regions_merger = Level_of_detail_facet_regions_merger_3<Kernel, Building>;
			
            Level_of_detail_coplanar_walls_merger_step_7(Buildings &buildings) :
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(1000))
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
            const FT m_tolerance;

            void process_building(Building &building) const {

                Walls &walls = building.walls;
                Output_facets output_facets;

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(output_facets);

                set_final_facets(output_facets, walls);
            }

            void set_final_facets(const Output_facets &output_facets, Walls &walls) const {

                walls.clear();
                walls.resize(output_facets.size());

                for (size_t i = 0; i < output_facets.size(); ++i)
                    set_final_facet(output_facets[i], walls[i]);
            }

            void set_final_facet(const Output_facet &output_facet, Wall &wall) const {

                Boundary &boundary = wall.boundary;
                boundary.clear();

                const size_t n = output_facet.size();
                for (size_t i = 0; i < n; ++i) {
                    
                    const size_t im = (i + n -1) % n;
                    const size_t ip = (i + 1) % n;

                    const Point_3 &p1 = output_facet[im];
                    const Point_3 &p2 = output_facet[i];
                    const Point_3 &p3 = output_facet[ip];

                    if (!are_collinear(p1, p2, p3)) boundary.push_back(p2);
                }
            }

            bool are_collinear(const Point_3 &p1, const Point_3 &p2, const Point_3 &p3) const {

                const Triangle_3 triangle = Triangle_3(p1, p2, p3);
                const FT area = CGAL::abs(static_cast<FT>(std::sqrt(CGAL::to_double(triangle.squared_area()))));

                return area < m_tolerance;
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COPLANAR_WALLS_MERGER_STEP_7_H