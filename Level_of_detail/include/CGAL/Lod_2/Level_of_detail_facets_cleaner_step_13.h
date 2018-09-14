#ifndef CGAL_LEVEL_OF_DETAIL_FACETS_CLEANER_STEP_13_H
#define CGAL_LEVEL_OF_DETAIL_FACETS_CLEANER_STEP_13_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_facets_cleaner_step_13 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;

            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron   = typename Building::Polyhedron;
			using Polyhedrons  = typename Building::Polyhedrons;
            
			using Clean_facet  = typename Building::Clean_facet;
			using Clean_facets = typename Building::Clean_facets;

			using Polyhedron_facet    = typename Polyhedron::Facet;
			using Polyhedron_facets   = typename Polyhedron::Facets;
			using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Output_regions = typename Building::Output_regions;

            Level_of_detail_facets_cleaner_step_13(Buildings &buildings) :
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(100000))
            { }

            void create_clean_facets() {

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
            const FT   m_tolerance;
            
            void process_building(Building &building) const {
            
                Clean_facets &clean_facets = building.clean_facets;
                clean_facets.clear();
                
                building.is_clean = true;
                const Polyhedrons &polyhedrons = building.polyhedrons;
                
				for (size_t i = 0; i < polyhedrons.size(); ++i)
					process_polyhedron(i, polyhedrons, clean_facets);

                Output_regions &output_regions = building.output_regions;
                output_regions.clear();

                output_regions.push_back(clean_facets);
            }

            void process_polyhedron(const size_t polyhedron_index, const Polyhedrons &polyhedrons, Clean_facets &clean_facets) const {

                const Polyhedron &polyhedron = polyhedrons[polyhedron_index];
                Clean_facet clean_facet;

                if (!polyhedron.is_valid) return;

				const Polyhedron_facets   &facets   = polyhedron.facets;
				const Polyhedron_vertices &vertices = polyhedron.vertices;

                for (size_t i = 0; i < facets.size(); ++i) {
                    const Polyhedron_facet &facet = facets[i];

                    if (!facet.is_valid) continue;

                    clean_facet.clear();
                    clean_facet.resize(facet.indices.size());

                    for (size_t j = 0; j < facet.indices.size(); ++j)
                        clean_facet[j] = vertices[facet.indices[j]];

                    if (!is_interior_facet(clean_facet, polyhedron_index, polyhedrons)) clean_facets.push_back(clean_facet);
                }
            }

            bool is_interior_facet(const Clean_facet &clean_facet, const size_t polyhedron_index, const Polyhedrons &polyhedrons) const {

                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    const Polyhedron &polyhedron = polyhedrons[i];
                    
                    if (!polyhedron.is_valid || i == polyhedron_index) continue;

                    const Polyhedron_facets   &facets   = polyhedron.facets;
				    const Polyhedron_vertices &vertices = polyhedron.vertices;

                    for (size_t j = 0; j < facets.size(); ++j) {
                        if (!facets[j].is_valid) continue;

                        if (are_equal_facets(clean_facet, facets[j], vertices)) 
                            return true;
                    }
                }
                return false;
            }

            bool are_equal_facets(const Clean_facet &clean_facet, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices) const {

                if (clean_facet.size() != facet.indices.size()) return false;

                size_t count = 0;
                for (size_t i = 0; i < clean_facet.size(); ++i) {
                    for (size_t j = 0; j < facet.indices.size(); ++j) {
                        
                        if (are_equal_points(clean_facet[i], vertices[facet.indices[j]])) {
                            
                            ++count;
                            break;
                        }
                    }
                }
                return count == clean_facet.size();
            }

            bool are_equal_points(const Point_3 &p, const Point_3 &q) const {
                return CGAL::abs(p.x() - q.x()) < m_tolerance && CGAL::abs(p.y() - q.y()) < m_tolerance && CGAL::abs(p.z() - q.z()) < m_tolerance;
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_FACETS_CLEANER_STEP_13_H