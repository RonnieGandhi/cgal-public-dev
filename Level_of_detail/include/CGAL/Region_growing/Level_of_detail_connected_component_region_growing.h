#ifndef CGAL_LEVEL_OF_DETAIL_CONNECTED_COMPONENT_REGION_GROWING_H
#define CGAL_LEVEL_OF_DETAIL_CONNECTED_COMPONENT_REGION_GROWING_H

// STL includes.
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData>
		class Level_of_detail_connected_component_region_growing {

		public:
			using Kernel = InputKernel;
            using Input  = InputData;

			using FT 	  = typename Kernel::FT;
			using Point_3 = typename Kernel::Point_3;
	
			using Index   = int;
			using Indices = std::vector<Index>;
			
            using Components = std::vector<Indices>;

			using Point_with_index = std::pair<Point_3, Index>;
            using Neighbours       = std::vector<Point_with_index>;

			using Point_map = CGAL::First_of_pair_property_map<Point_with_index>;

            using Search_traits_3 = CGAL::Search_traits_3<Kernel>;
            using Search_traits   = CGAL::Search_traits_adapter<Point_with_index, Point_map, Search_traits_3>;
			using Fuzzy_tree      = CGAL::Kd_tree<Search_traits>;
            using Fuzzy_sphere    = CGAL::Fuzzy_sphere<Search_traits>;

            using States = std::map<Index, bool>;
            using Queue  = std::queue<Index>;

			Level_of_detail_connected_component_region_growing(const Input &input, const Indices &indices) : 
            m_input(input),
            m_indices(indices),
            m_scale(FT(1)),
            m_min_size(10)
			{ }

            void set_scale(const FT new_value) {
                m_scale = new_value;
            }

            void set_min_points(const size_t new_value) {
                m_min_size = new_value;
            }

			void apply(Components &components) const {				
				
                Fuzzy_tree tree;
                components.clear();

				create_tree(tree);
                get_connected_components(tree, components);
			}

		private:
            const Input   &m_input;
            const Indices &m_indices;

            FT     m_scale;
            size_t m_min_size;

            void create_tree(Fuzzy_tree &tree) const {

                for (size_t i = 0; i < m_indices.size(); ++i)
                    tree.insert(std::make_pair(m_input.point(m_indices[i]), m_indices[i]));
            }

            void get_connected_components(const Fuzzy_tree &tree, Components &components) const {

				Queue   queue;
                Indices component;
                States  states;
                
                for (size_t i = 0; i < m_indices.size(); ++i)
                    states[m_indices[i]] = false;

				for (size_t i = 0; i < m_indices.size(); ++i) {
                    const Index point_index = m_indices[i];
					
					if (states.at(point_index)) 
                        continue;

					component.clear();
					component.push_back(point_index);
					
					states[point_index] = true;

					for (size_t j = 0; j < queue.size(); ++j) queue.pop();
					queue.push(point_index);

					propagate(tree, queue, component, states);
					if (component.size() > m_min_size) components.push_back(component);
				}
            }

            void propagate(const Fuzzy_tree &tree, Queue &queue, Indices &component, States &states) const {

				while (!queue.empty()) {
					
					const Index point_index = queue.front();
					queue.pop();

                    const Fuzzy_sphere sphere = Fuzzy_sphere(m_input.point(point_index), m_scale);
                    
                    Neighbours neighbours;
                    find_nearest_neighbours(tree, sphere, neighbours);
					
					for (size_t i = 0; i < neighbours.size(); ++i) {
						const Index neighbour_index = neighbours[i].second;

						if (!states.at(neighbour_index)) {

							queue.push(neighbour_index);
							component.push_back(neighbour_index);
							states[neighbour_index] = true;
						}
					}
				}
			}

            void find_nearest_neighbours(const Fuzzy_tree &tree, const Fuzzy_sphere &sphere, Neighbours &neighbours) const {

				neighbours.clear();
				tree.search(std::back_inserter(neighbours), sphere);
			}
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_CONNECTED_COMPONENT_REGION_GROWING_H