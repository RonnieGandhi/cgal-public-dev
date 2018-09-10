#ifndef CGAL_LEVEL_OF_DETAIL_WEIGHT_QUALITY_POLYHEDRON_ESTIMATOR_STEP_10_H
#define CGAL_LEVEL_OF_DETAIL_WEIGHT_QUALITY_POLYHEDRON_ESTIMATOR_STEP_10_H

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_weight_quality_polyhedron_estimator_step_10 {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;
            
            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Index   = typename Building::Index;
            using Indices = typename Building::Indices;

            using Polyhedron_facet    = typename Polyhedron::Facet;
            using Polyhedron_facets   = typename Polyhedron::Facets;
            using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Graphcut_facet  = typename Building::Graphcut_facet;
            using Graphcut_facets = typename Building::Graphcut_facets;

            using Graphcut_facet_data      = typename Graphcut_facet::Data;
            using Graphcut_facet_data_pair = typename Graphcut_facet::Data_pair;

            Level_of_detail_weight_quality_polyhedron_estimator_step_10(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height),
            m_tolerance(FT(1) / FT(100000))
            { }

            void estimate() {

                if (m_buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid && building.interior_indices.size() != 0)
                        process_building(building);
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;

            const FT m_ground_height;
            const FT m_tolerance;

            void process_building(Building &building) {

                Graphcut_facets &gc_facets = building.graphcut_facets;
                create_graphcut_facets(building, gc_facets);
            }

            void create_graphcut_facets(const Building &building, Graphcut_facets &gc_facets) const {
                
                gc_facets.clear();
                
                const Indices &interior_indices = building.interior_indices;
                const Polyhedrons &polyhedrons  = building.polyhedrons;

                for (int i = 0; i < polyhedrons.size(); ++i)
                    add_graphcut_facets(interior_indices, polyhedrons, i, gc_facets);
            }

            void add_graphcut_facets(const Indices &interior_indices, const Polyhedrons &polyhedrons, const int poly_index, Graphcut_facets &gc_facets) const {

                for (int i = 0; i < polyhedrons[poly_index].facets.size(); ++i)
                    process_polyhedron_facet(interior_indices, polyhedrons, poly_index, i, gc_facets);
            }

            void process_polyhedron_facet(const Indices &interior_indices, const Polyhedrons &polyhedrons, const int poly_index, const int facet_index, Graphcut_facets &gc_facets) const {

                if (was_already_added(poly_index, facet_index, gc_facets))
                    return;

                const Polyhedron &polyhedron = polyhedrons[poly_index];
                
                const Polyhedron_facets   &facets   = polyhedron.facets;
                const Polyhedron_vertices &vertices = polyhedron.vertices;

                Graphcut_facet gc_facet;
                gc_facet.data = std::make_pair(poly_index, facet_index);

                gc_facet.weight  = compute_weight(facets[facet_index], vertices);
                gc_facet.quality = compute_quality(interior_indices, facets[facet_index], vertices);

                find_neighbours(polyhedrons, facets, vertices, poly_index, facet_index, gc_facet.neighbours);
                gc_facets.push_back(gc_facet);
            }

            bool was_already_added(const int poly_index, const int facet_index, const Graphcut_facets &gc_facets) const {

                for (size_t i = 0; i < gc_facets.size(); ++i) {
                    const Graphcut_facet_data_pair &data_pair = gc_facets[i].neighbours;

                    const Graphcut_facet_data &neigh_1 = data_pair.first;
                    const Graphcut_facet_data &neigh_2 = data_pair.second;

                    if (neigh_1.first == poly_index && neigh_1.second == facet_index) return true;
                    if (neigh_2.first == poly_index && neigh_2.second == facet_index) return true;
                }
                return false;
            }

            FT compute_weight(const Polyhedron_facet &, const Polyhedron_vertices &) const {
                return FT(1);
            }

            FT compute_quality(const Indices &interior_indices, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices) const {
                return FT(0);
            }

            void find_neighbours(const Polyhedrons &polyhedrons, const Polyhedron_facets &facets, const Polyhedron_vertices &vertices, 
            const int poly_index, const int facet_index, Graphcut_facet_data_pair &neighbours) const {

                Graphcut_facet_data &neigh_1 = neighbours.first;
                Graphcut_facet_data &neigh_2 = neighbours.second;

                neigh_1.first  = poly_index;
                neigh_1.second = facet_index;

                find_neighbour(polyhedrons, poly_index, facets[facet_index], vertices, neigh_2);
            }

            void find_neighbour(const Polyhedrons &polyhedrons, const int poly_index, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, Graphcut_facet_data &neigh) const {

                neigh.first  = -1;
                neigh.second = -1;

                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    if (i == poly_index) continue;

                    for (size_t j = 0; j < polyhedrons[i].facets.size(); ++j) {
                        if (are_equal_facets(facet, vertices, polyhedrons[i].facets[j], polyhedrons[i].vertices)) {

                            neigh.first  = i;
                            neigh.second = j;

                            return;
                        }
                    }
                }
            }

            bool are_equal_facets(const Polyhedron_facet &f1, const Polyhedron_vertices &v1, const Polyhedron_facet &f2, const Polyhedron_vertices &v2) const {

                if (f1.indices.size() != f2.indices.size()) return false;

                size_t count = 0;
                for (size_t i = 0; i < f1.indices.size(); ++i) {
                    for (size_t j = 0; j < f2.indices.size(); ++j) {

                        if (CGAL::abs(v1[f1.indices[i]].x() - v2[f2.indices[j]].x()) < m_tolerance && 
                            CGAL::abs(v1[f1.indices[i]].y() - v2[f2.indices[j]].y()) < m_tolerance && 
                            CGAL::abs(v1[f1.indices[i]].z() - v2[f2.indices[j]].z()) < m_tolerance) {
                            ++count; break;
                        }
                    }
                }
                return count == f1.indices.size();
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_WEIGHT_QUALITY_POLYHEDRON_ESTIMATOR_STEP_10_H