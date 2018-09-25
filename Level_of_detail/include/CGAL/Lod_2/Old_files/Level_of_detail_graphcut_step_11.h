#ifndef CGAL_LEVEL_OF_DETAIL_GRAPHCUT_STEP_11_H
#define CGAL_LEVEL_OF_DETAIL_GRAPHCUT_STEP_11_H 

// STL includes.
#include <map>
#include <vector>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace Maxflow {
	#include <CGAL/internal/auxiliary/graph.h>
}

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_graphcut_step_11 {

		public:
            using Kernel    = InputKernel;
			using Building  = InputBuilding;
			using Buildings = InputBuildings;

			using FT = typename Kernel::FT;

			using Graph   = Maxflow::Graph;
			using Node_id = typename Graph::node_id;

			using Buildings_iterator = typename Buildings::iterator;

			using Polyhedron  = typename Building::Polyhedron;
			using Polyhedrons = typename Building::Polyhedrons;
			
			using Graphcut_facet  = typename Building::Graphcut_facet;
			using Graphcut_facets = typename Building::Graphcut_facets;

			using Neighbour  = typename Graphcut_facet::Data;
			using Neighbours = typename Graphcut_facet::Data_pair;

			Level_of_detail_graphcut_step_11() : 
			m_beta(FT(0.1)),
			m_run_test(false)
			{ }

			void set_beta(const FT beta) {
				
			}

			void solve(Buildings &buildings) const {

				if (m_run_test) 
					run_test();

				if (buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    Building &building = bit->second;

					if (building.is_valid && building.polyhedrons.size() != 0) process_building(building);
					else building.is_valid = false;
                }
			}

		private:
			FT m_beta;
			const bool m_run_test;

			void run_test() const {

				Graph::node_id nodes[2];
				Graph *g = new Graph();

				nodes[0] = g -> add_node();
				nodes[1] = g -> add_node();

				g->set_tweights(nodes[0], 1, 5);
				g->set_tweights(nodes[1], 2, 6);

				g->add_edge(nodes[0], nodes[1], 3, 4);

				Graph::flowtype flow = g->maxflow();

				printf("Flow = %d\n", flow);
				printf("Minimum cut:\n");
				if (g->what_segment(nodes[0]) == Graph::SOURCE)
					printf("node0 is in the SOURCE set\n");
				else
					printf("node0 is in the SINK set\n");
				if (g->what_segment(nodes[1]) == Graph::SOURCE)
					printf("node1 is in the SOURCE set\n");
				else
					printf("node1 is in the SINK set\n");

				delete g;
				exit(1);
			}

			void process_building(Building &building) const {

				Polyhedrons 		  &polyhedrons = building.polyhedrons;
				Graphcut_facets &graphcut_facets   = building.graphcut_facets;

				const unsigned int num_nodes = polyhedrons.size();

				Node_id *pNodes = new Node_id[num_nodes];
				Graph    *graph = new Graph();

				set_graph_nodes(polyhedrons, pNodes, graph);
				set_graph_edges(graphcut_facets, pNodes, graph);

				graph->maxflow();
				set_solution(pNodes, graph, polyhedrons);

				delete   graph;
				delete[] pNodes;
            }

			void set_graph_nodes(Polyhedrons &polyhedrons, Node_id *pNodes, Graph *graph) const {

				/*
				FT total_weight = FT(0);
				for (size_t i = 0; i < polyhedrons.size(); ++i) {
					
					Polyhedron &polyhedron = polyhedrons[i];
					total_weight += polyhedron.weight;
				}

				for (size_t i = 0; i < polyhedrons.size(); ++i)
					polyhedrons[i].weight /= total_weight; */

				for (size_t i = 0; i < polyhedrons.size(); ++i) {
					pNodes[i] = graph->add_node();

					const Polyhedron &polyhedron = polyhedrons[i];

					const FT in  = polyhedron.in;
					const FT out = polyhedron.out;

					// std::cout << "nodes: " << in << " : " << out << std::endl;

					CGAL_precondition(in  >= FT(0) && in  <= FT(1));
					CGAL_precondition(out >= FT(0) && out <= FT(1));

					const FT node_weight = polyhedron.weight;
					CGAL_precondition(node_weight >= FT(0));

					const FT cost_in  = get_graph_node_cost(in , node_weight);
					const FT cost_out = get_graph_node_cost(out, node_weight);

					graph->add_tweights(pNodes[i], CGAL::to_double(cost_in), CGAL::to_double(cost_out));
				}
			}

			FT get_graph_node_cost(const FT node_value, const FT node_weight) const {
				return node_weight * node_value;
			} 

			void set_graph_edges(Graphcut_facets &graphcut_facets, const Node_id *pNodes, Graph *graph) const {

				/*
				FT total_weight = FT(0);
				for (size_t i = 0; i < graphcut_facets.size(); ++i) {
					
					Graphcut_facet &graphcut_facet = graphcut_facets[i];
					total_weight += graphcut_facet.weight;
				}

				for (size_t i = 0; i < graphcut_facets.size(); ++i)
					graphcut_facets[i].weight /= total_weight; */

				for (size_t i = 0; i < graphcut_facets.size(); ++i) {
					
					const Graphcut_facet &graphcut_facet = graphcut_facets[i];
					const Neighbours 	 &neighbours 	 = graphcut_facet.neighbours;

					const Neighbour &neigh_1 = neighbours.first;
					const Neighbour &neigh_2 = neighbours.second;

					const int polyhedron_index_1 = neigh_1.first;
					const int polyhedron_index_2 = neigh_2.first;

					if (polyhedron_index_1 < 0 || polyhedron_index_2 < 0) continue;

					const FT edge_weight  = graphcut_facet.weight;
					const FT edge_quality = graphcut_facet.quality;

					CGAL_precondition(edge_weight  >= FT(0));
					CGAL_precondition(edge_quality >= FT(0) && edge_quality <= FT(1));

					add_graph_edge(pNodes, polyhedron_index_1, polyhedron_index_2, edge_weight, edge_quality, graph);
				}
			}

			void add_graph_edge(
				const Node_id *pNodes, 
				const int i,
				const int j,
				const FT edge_weight, 
				const FT edge_quality, 
				Graph *graph) const {

				const FT cost_value = get_graph_edge_cost(edge_weight, edge_quality);
				graph->add_edge(pNodes[i], pNodes[j], CGAL::to_double(cost_value), CGAL::to_double(cost_value));
			}

			FT get_graph_edge_cost(const FT edge_weight, const FT edge_quality) const {
				// std::cout << edge_weight << std::endl;
				return m_beta * edge_weight; // * inverse(edge_quality);
			}

			FT inverse(const FT value) const {
				
				CGAL_precondition(value >= FT(0) && value <= FT(1));
				return FT(1) - value;
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

#endif // CGAL_LEVEL_OF_DETAIL_GRAPHCUT_STEP_11_H