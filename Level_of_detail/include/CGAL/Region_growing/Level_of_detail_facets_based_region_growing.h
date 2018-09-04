#ifndef CGAL_LEVEL_OF_DETAIL_FACETS_BASED_REGION_GROWING_H
#define CGAL_LEVEL_OF_DETAIL_FACETS_BASED_REGION_GROWING_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif

// STL includes.
#include <map>
#include <queue>
#include <cmath>
#include <time.h>
#include <vector>
#include <stdio.h>
#include <cassert>
#include <stdlib.h>
#include <algorithm>

// CGAL includes.
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel>
		class Level_of_detail_facets_based_region_growing {

		public:
			using Kernel = InputKernel;

			using FT 	   = typename Kernel::FT;
			using Point_3  = typename Kernel::Point_3;
			using Vector_3 = typename Kernel::Vector_3;

			using Color = CGAL::Color;

			using Region_facet   = std::vector<Point_3>;
            using Region_facets  = std::vector<Region_facet>;
			using Output_regions = std::vector<Region_facets>;
			using Output_colors  = std::vector<Color>;
			
			using Input_facet  = std::pair<Region_facet, Color>;
			using Input_facets = std::vector<Input_facet>;

			Input_facets input_facets;

			using Log = CGAL::LOD::Mylog;

			using Region  = std::vector<size_t>;
			using Regions = std::vector<Region>;

			using States = std::map<size_t, bool>;
			
			using Neighbours   = std::vector<size_t>;
			using Connectivity = std::map<size_t, Neighbours>;

			using Queue = std::queue<size_t>;

			Level_of_detail_facets_based_region_growing(const Input_facets &input_facets) : 
			m_input_facets(input_facets),
			m_no_rg(false),
			m_tolerance(FT(1) / FT(100000)) { 

				srand(time(NULL));
			}

			void create_regions(Output_regions &output_regions) const {
				
				output_regions.clear();
				if (m_input_facets.size() == 0) return;

				Regions regions;
				grow_regions(regions);

				Output_colors output_colors;
				create_output_regions(regions, output_regions, output_colors);

				print_output_regions(output_regions, output_colors);
			}

		private:
			const Input_facets &m_input_facets;
			const bool m_no_rg;

			const FT m_tolerance;

			void grow_regions(Regions &regions) const {

				if (m_no_rg) do_not_grow_all_regions(regions);
				else grow_all_regions(regions);
			}

			void do_not_grow_all_regions(Regions &regions) const {
				
				regions.clear();
				regions.resize(m_input_facets.size());

				for (size_t i = 0; i < regions.size(); ++i) {
					Region &region = regions[i];

					region.resize(1);
					region[0] = i;
				}
			}

			void create_output_regions(const Regions &regions, Output_regions &output_regions, Output_colors &output_colors) const {
				
				output_regions.clear();
				output_regions.resize(regions.size());

				output_colors.clear();
				output_colors.resize(regions.size());

				for (size_t i = 0; i < regions.size(); ++i) {
					
					output_regions[i].resize(regions[i].size());
					output_colors[i] = generate_random_color();

					for (size_t j = 0; j < regions[i].size(); ++j) {
						
						const size_t index = regions[i][j];
						create_output_facet(index, output_regions[i][j]);
					}
				}
			}

			void create_output_facet(const size_t facet_index, Region_facet &region_facet) const {

				region_facet.clear();
				region_facet.resize(m_input_facets[facet_index].first.size());

				for (size_t i = 0; i < m_input_facets[facet_index].first.size(); ++i)
					region_facet[i] = m_input_facets[facet_index].first[i];
			}

			void print_output_regions(const Output_regions &output_regions, const Output_colors &output_colors) const {

				// Log log;
				// log.print_rg_output_facets(output_regions, output_colors, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "6_region_growing");
			}

            Color generate_random_color() const {

				const int r = rand() % 255;
				const int g = rand() % 255;
				const int b = rand() % 255;

				return Color(r, g, b);
			}

			void grow_all_regions(Regions &regions) const {
				
				Connectivity connectivity;
				create_connectivity(connectivity);

				States states;
				create_default_states(states);

				apply_region_growing_algorithm(connectivity, states, regions);
			}

			void create_default_states(States &states) const {
				
				states.clear();
				for (size_t i = 0; i < m_input_facets.size(); ++i)
					states[i] = false;
			}

			void create_connectivity(Connectivity &connectivity) const {
				
				connectivity.clear();
				Neighbours neighbours;

				for (size_t i = 0; i < m_input_facets.size(); ++i) {
					
					find_facet_neighbours(m_input_facets[i], i, neighbours);
					connectivity[i] = neighbours;

					// std::cout << m_input_facets[i].first.size() << " : " << neighbours.size() << std::endl;
				}
			}

			void find_facet_neighbours(const Input_facet &input_facet, const size_t facet_index, Neighbours &neighbours) const {
				neighbours.clear();

				for (size_t i = 0; i < m_input_facets.size(); ++i)
					if (i != facet_index && share_an_edge(input_facet, m_input_facets[i]))
						neighbours.push_back(i);
			}

			bool share_an_edge(const Input_facet &f1, const Input_facet &f2) const {

				for (size_t i = 0; i < f1.first.size(); ++i) {
					const size_t ip = (i + 1) % f1.first.size();

					for (size_t j = 0; j < f2.first.size(); ++j) {
						const size_t jp = (j + 1) % f2.first.size();

						if (are_equal_edges(f1.first[i], f1.first[ip], f2.first[j], f2.first[jp])) 
							return true;

						// if (are_equal_points(f1.first[i], f2.first[j]))
						//     return true;
					}
				}
				return false;
			}

			bool are_equal_edges(const Point_3 &p1, const Point_3 &p2, const Point_3 &q1, const Point_3 &q2) const {
                return (are_equal_points(p1, q1) && are_equal_points(p2, q2)) || (are_equal_points(p1, q2) && are_equal_points(p2, q1));
            }

            bool are_equal_points(const Point_3 &p, const Point_3 &q) const {

                const FT eps = m_tolerance;
                return (CGAL::abs(p.x() - q.x()) < eps) && (CGAL::abs(p.y() - q.y()) < eps) && (CGAL::abs(p.z() - q.z()) < eps);
            }

			void apply_region_growing_algorithm(const Connectivity &connectivity, States &states, Regions &regions) const {
				regions.clear();

				Region region; Queue queue;
				for (size_t i = 0; i < m_input_facets.size(); ++i) {
					
					const Input_facet &input_facet = m_input_facets[i];
					if (states.at(i)) continue;

					region.clear();
					region.push_back(i);
					
					states[i] = true;

					for (size_t j = 0; j < queue.size(); ++j) queue.pop();
					queue.push(i);

					grow_region(connectivity, queue, region, states);
					regions.push_back(region);
				}
			}

			void grow_region(const Connectivity &connectivity, Queue &queue, Region &region, States &states) const {

				while (!queue.empty()) {
					
					const size_t facet_index = queue.front();
					queue.pop();

					const Neighbours &neighbours = connectivity.at(facet_index);
					for (size_t i = 0; i < neighbours.size(); ++i) {

						const size_t neighbour_index = neighbours[i];
						/* if (!states.at(neighbour_index) && m_input_facets[neighbour_index].first.size() == 3) {
							
							std::cout << m_input_facets[facet_index].first.size() << " : " << neighbours.size() << std::endl;
							std::cout << "cond: " << satisfies_local_conditions(facet_index, neighbour_index) << std::endl;
						} */

						if (!states.at(neighbour_index) && satisfies_local_conditions(facet_index, neighbour_index)) {

							queue.push(neighbour_index);
							region.push_back(neighbour_index);
							states[neighbour_index] = true;
						}
					}
				}
			}

			bool satisfies_local_conditions(const size_t f1_index, const size_t f2_index) const {	
				return /* are_within_angle_tolerance(f1_index, f2_index) && */ are_coplanar(f1_index, f2_index);
			}

			/*
			bool are_within_angle_tolerance(const size_t f1_index, const size_t f2_index) const {
				
				Vector_3 n1, n2;
				get_normal(f1_index, n1);
				get_normal(f2_index, n2);

				const FT scalar = CGAL::scalar_product(n1, n2);
				return (FT(1) - CGAL::abs(scalar)) < m_tolerance;
			}

			Vector_3 get_normal(const size_t facet_index, Vector_3 &normal) const {

				const Input_facet &input_facet = m_input_facets[facet_index];
				CGAL_precondition(input_facet.first.size() >= 3);

				const Point_3 &p1 = input_facet.first[0];
				const Point_3 &p2 = input_facet.first[1];
				const Point_3 &p3 = input_facet.first[2];

				const Vector_3 v1 = Vector_3(p1, p2);
				const Vector_3 v2 = Vector_3(p1, p3);

				const Vector_3 cross = CGAL::cross_product(v1, v2);
				normal = cross / static_cast<FT>(CGAL::sqrt(CGAL::to_double(cross.squared_length())));
			} */

			bool are_coplanar(const size_t f1_index, const size_t f2_index) const {

				const Input_facet &v1 = m_input_facets[f1_index];
				const Input_facet &v2 = m_input_facets[f2_index];

				size_t count = 0;
				for (size_t i = 0; i < v1.first.size(); ++i) {

					const size_t ip  = (i + 1) % v1.first.size();
					const size_t ipp = (i + 2) % v1.first.size();

					for (size_t j = 0; j < v2.first.size(); ++j)
						if (is_coplanar(v1.first[i], v1.first[ip], v1.first[ipp], v2.first[j])) ++count;
				}
				return count == v1.first.size() * v2.first.size();
			}

			bool is_coplanar(const Point_3 &p1, const Point_3 &p2, const Point_3 &p3, const Point_3 &p4) const {

				const Vector_3 v1 = Vector_3(p1, p2);
				const Vector_3 v2 = Vector_3(p1, p3);
				const Vector_3 v3 = Vector_3(p1, p4);

				const Vector_3 v4 =  CGAL::cross_product(v2, v3);
				const FT result   = CGAL::scalar_product(v1, v4);

				return CGAL::abs(result) < m_tolerance;
			}
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_FACETS_BASED_REGION_GROWING_H