#ifndef CGAL_LEVEL_OF_DETAIL_GRAPHCUT_3_H
#define CGAL_LEVEL_OF_DETAIL_GRAPHCUT_3_H 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace Maxflow {
	#include <CGAL/internal/auxiliary/graph.h>
}

// New CGAL includes.
#include <CGAL/Graphcut/Graphcut_3/Level_of_detail_polyhedron_in_out_estimator.h>
#include <CGAL/Graphcut/Graphcut_3/Level_of_detail_facets_weight_quality_estimator.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding, class InputBuildings>
		class Level_of_detail_graphcut_3 {

		public:
            using Kernel    = InputKernel;
			using Input     = InputData;
			using Building  = InputBuilding;
			using Buildings = InputBuildings;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_3 Point_3;

			typedef Maxflow::Graph 			Graph;
			typedef typename Graph::node_id Node_id;

			using Buildings_iterator = typename Buildings::iterator;

			using Polyhedron  = typename Building::Polyhedron;
			using Polyhedrons = typename Building::Polyhedrons;
			
			using Graphcut_facet  = typename Building::Graphcut_facet;
			using Graphcut_facets = typename Building::Graphcut_facets;

			using Polyhedron_in_out_estimator 	  = CGAL::LOD::Level_of_detail_polyhedron_in_out_estimator<Kernel, Input, Building>;
			using Facets_weight_quality_estimator = CGAL::LOD::Level_of_detail_facets_weight_quality_estimator<Kernel, Input, Building>;

			Level_of_detail_graphcut_3(const Input &input, const FT ground_height) : 
			m_input(input),
			m_ground_height(ground_height),
			m_alpha(FT(1)),
			m_beta(FT(100000)),
			m_gamma(FT(1000)) { 

				
			}

			void set_alpha(const FT alpha) {
				// m_alpha = alpha;
			}

			void set_beta(const FT beta) {
				// m_beta = beta;
			}

			void set_gamma(const FT gamma) {
				// m_gamma = gamma;
			}

			void solve(Buildings &buildings) const {

				if (buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    Building &building = bit->second;

					if (building.is_valid && building.shapes.size() != 0 && building.polyhedrons.size() != 0) 
						process_building(building);
					else 
						building.is_valid = false;
                }
			}

		private:
			const Input &m_input;
			const FT m_ground_height;

			FT m_alpha;
			FT m_beta;
			FT m_gamma;

			void process_building(Building &building) const {

				Polyhedron_in_out_estimator polyhedron_in_out_estimator(m_input, m_ground_height);
				polyhedron_in_out_estimator.estimate(building);

				Facets_weight_quality_estimator facets_weight_quality_estimator(m_input, m_ground_height);
				facets_weight_quality_estimator.estimate(building);

				Polyhedrons &polyhedrons 	 	 	   = building.polyhedrons;
				const Graphcut_facets &graphcut_facets = building.graphcut_facets;

				const unsigned int num_nodes = polyhedrons.size();

				Node_id *pNodes = new Node_id[num_nodes];
				Graph    *graph = new Graph();

				set_graph_nodes(polyhedrons, pNodes, graph);
				set_graph_edges(graphcut_facets, pNodes, graph);

				// graph->maxflow();
				
				set_polyhedrons_in_out(polyhedrons);

				// set_solution(pNodes, graph, polyhedrons);

				delete   graph;
				delete[] pNodes;
            }

			void set_graph_nodes(const Polyhedrons &polyhedrons, Node_id *pNodes, Graph *graph) const {

				for (size_t i = 0; i < polyhedrons.size(); ++i) {
					pNodes[i] = graph->add_node();

					const Polyhedron &polyhedron = polyhedrons[i];

					const FT in  = polyhedron.in;
					const FT out = polyhedron.out;

					CGAL_precondition(in  >= FT(0) && in  <= FT(1));
					CGAL_precondition(out >= FT(0) && out <= FT(1));

					const FT cost_in  = get_graph_node_cost(in, polyhedron);
					const FT cost_out = get_graph_node_cost(out, polyhedron);

					graph->add_tweights(pNodes[i], CGAL::to_double(cost_in), CGAL::to_double(cost_out));
				}
			}

			FT get_graph_node_cost(const FT node_value, const Polyhedron &polyhedron) const {
				
				return get_graph_node_weight(polyhedron) * node_value;
			} 

			FT get_graph_node_weight(const Polyhedron &polyhedron) const {
				
				return m_beta * get_polyhedron_weight(polyhedron);
			}

			FT get_polyhedron_weight(const Polyhedron &polyhedron) const {
				
				return FT(1); // compute volume of the polyhedron here!
			}

			void set_graph_edges(const Graphcut_facets &graphcut_facets, const Node_id *pNodes, Graph *graph) const {

				for (size_t i = 0; i < graphcut_facets.size(); ++i) {
					const Graphcut_facet &graphcut_facet = graphcut_facets[i];

					const size_t polyhedron_index_1 = graphcut_facet.neighbours.first;
					const size_t polyhedron_index_2 = graphcut_facet.neighbours.second;

					const FT edge_weight  = graphcut_facet.weight;
					const FT edge_quality = graphcut_facet.quality;

					CGAL_precondition(edge_weight  >= FT(0));
					CGAL_precondition(edge_quality >= FT(0));

					add_graph_edge(pNodes, polyhedron_index_1, polyhedron_index_2, edge_weight, edge_quality, graph);
				}
			}

			void add_graph_edge(
				const Node_id *pNodes, 
				const size_t i,
				const size_t j,
				const FT edge_weight, 
				const FT edge_quality, 
				Graph *graph) const {

				const FT cost_value = get_graph_edge_cost(edge_weight, edge_quality);
				graph->add_edge(pNodes[i], pNodes[j], CGAL::to_double(cost_value), CGAL::to_double(cost_value));
			}

			FT get_graph_edge_cost(const FT edge_weight, const FT edge_quality) const {
				return edge_weight * edge_quality;
			}

			void set_polyhedrons_in_out(Polyhedrons &polyhedrons) const {

				for (size_t i = 0; i < polyhedrons.size(); ++i) {
					Polyhedron &polyhedron = polyhedrons[i];

					if (polyhedron.in > FT(1) / FT(2)) { polyhedron.is_valid = true;  continue; }
					else { polyhedron.is_valid = false; continue; }
				}
			}

			void set_solution(const Node_id *pNodes, Graph *graph, Polyhedrons &polyhedrons) const {

				for (size_t i = 0; i < polyhedrons.size(); ++i) {
					Polyhedron &polyhedron = polyhedrons[i];

					if (graph->what_segment(pNodes[i]) == Graph::SOURCE) { polyhedron.is_valid = true;  continue; }
					if (graph->what_segment(pNodes[i]) == Graph::SINK)   { polyhedron.is_valid = false; continue; }

					std::cout << "Error: graphcut 3: cannot be here, smth is wrong!" << std::endl;
					exit(0);
				}
			}
		};

	} // CGAL

} // LOD

#endif // CGAL_LEVEL_OF_DETAIL_GRAPHCUT_3_H