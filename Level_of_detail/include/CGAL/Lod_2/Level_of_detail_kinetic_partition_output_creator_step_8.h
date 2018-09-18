#ifndef CGAL_LEVEL_OF_DETAIL_KINETIC_PARTITION_OUTPUT_CREATOR_STEP_8_H
#define CGAL_LEVEL_OF_DETAIL_KINETIC_PARTITION_OUTPUT_CREATOR_STEP_8_H

// STL includes.
#include <set>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

// Jean Philippe includes.
#include <CGAL/Buildings/jean_philippe/propagation.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_kinetic_partition_output_creator_step_8 {
            
        public:
            using Kernel    = InputKernel;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using Buildings_iterator = typename Buildings::iterator;

            using JP_polygon  = typename Building::JP_polygon;
            using JP_polygons = typename Building::JP_polygons;

            using JP_kinetic_propagation = JPTD::Kinetic_Propagation;

            using JP_polyhedron  = JPTD::Partition_Polyhedron;
            using JP_polyhedrons = std::list<JP_polyhedron*>;

            using JP_polyhedrons_iterator = typename JP_polyhedrons::const_iterator;

            using JP_vertex = JPTD::Partition_Vertex;

            using JP_sequence  = std::list<JP_vertex*>;
            using JP_sequences = std::list<JP_sequence>;

            using JP_sequences_iterator = typename JP_sequences::const_iterator;
            using JP_sequence_set = std::set<JP_vertex*>;

            using JP_sequence_iterator     = typename JP_sequence::const_iterator;
            using JP_sequence_set_iterator = typename JP_sequence_set::const_iterator;

            using JP_conversions = std::map<const JP_vertex*, int>;

            using JP_side  = JPTD::Partition_Side;
            using JP_sides = std::list<JP_side>;

            using JP_sides_iterator = typename JP_sides::const_iterator;

            using JP_facet          = JPTD::Partition_Facet;
            using JP_facet_vertices = std::vector<JP_vertex*>;
            
            using JP_facets = std::list<JP_facet*>;
            using JP_facets_iterator = typename JP_facets::const_iterator;

            using JP_point_3 = JPTD::CGAL_Point_3;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Polyhedron_vertex   = typename Polyhedron::Vertex;
            using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Polyhedron_facet  = typename Polyhedron::Facet;
            using Polyhedron_facets = typename Polyhedron::Facets;
            
            using Log = CGAL::LOD::Mylog;

            typename Kernel::Compute_squared_distance_3 squared_distance_3;

            Level_of_detail_kinetic_partition_output_creator_step_8(Buildings &buildings) :
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(1000))
            { }

            void create() const {
                
                // debug_buildings();

                if (m_buildings.size() == 0) 
                    return; 
                    
                int count = 0; std::cout << std::endl; bool is_valid_building = true;
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it, ++count) {
                    
                    Building &building = b_it->second; is_valid_building = building.is_valid;
                    if (building.jp_polygons.size() == 0) is_valid_building = false;

                    building.index = count;
                    // std::cout << "index: " << building.index << " ";

                    std::cout << "percents: " << FT(count) / FT(m_buildings.size()) * FT(100) << "%" << std::endl;

                    // Log log; 
                    // log.save_only_convex_polygons(building.jp_polygons, "tmp/buildings/debug_building_" + std::to_string(building.index));

					if ( is_valid_building /* && 
                    building.index != 112 && 
                    building.index != 113 && 
                    building.index != 114 && 
                    building.index != 147 && 
                    building.index != 148 && 
                    building.index != 374 && 
                    building.index != 404 && 
                    building.index != 392*/ ) process_building(building);
                }
                std::cout << std::endl;
            }

        private:
            Buildings &m_buildings;
            const FT m_tolerance;

            void process_building(Building &building) const {

                const JP_polygons &jp_polygons = building.jp_polygons;
                JP_kinetic_propagation kinetic(jp_polygons);

	            if (!kinetic.data()) return;
                kinetic.run();

                set_output(kinetic, building);
                kinetic.delete_kinetic_data_structure();
            }

            void debug_buildings() const {
                
                // bad buildings are 82, 117, 148, 375

                // debug_building(10);
                
                // debug_building(82);
                // frame #1: 0x00000001005cd2da lod`JPTD::Partition::remove_bivalent_vertices(this=0x0000000117e5f380) at partition.cpp:775
                // 775                                  V.erase(V_map[v->id]);

                // debug_building(117);
                // frame #1: 0x00000001005cd2da lod`JPTD::Partition::remove_bivalent_vertices(this=0x0000000117e5f380) at partition.cpp:775
                // 775                                  V.erase(V_map[v->id]);

                // debug_building(148);
                // frame #1: 0x00000001005cd2da lod`JPTD::Partition::remove_bivalent_vertices(this=0x0000000117e5f380) at partition.cpp:775
                // 775                                  V.erase(V_map[v->id]);

                // debug_building(375);
                // frame #1: 0x00000001005cd2da lod`JPTD::Partition::remove_bivalent_vertices(this=0x0000000117e5f380) at partition.cpp:775
                // 775                                  V.erase(V_map[v->id]);

                // exit(1);
            }

            void debug_building(const size_t building_index) const {

                const std::string path = "/Users/danisimo/Documents/pipeline/logs/tmp/bad_input/debug_building_" + std::to_string(building_index) + ".ply";
                std::cout << "input file: " << path << std::endl;
 
                const std::string str1 = "./kinetic";
                const std::string str2 = "--input";
                const std::string str3 = "--check";
                const std::string str4 = "--print-drawings";

                char *val1 = strdup(str1.c_str());
                char *val2 = strdup(str2.c_str());
                char *val3 = strdup(path.c_str());
                char *val4 = strdup(str3.c_str());
                char *val5 = strdup(str4.c_str());

                int argc = 5;
                char *n_argv[] = { val1, val2, val3, val4, val5 };
                char **argv = n_argv;

                JP_kinetic_propagation kinetic(argc, argv);

	            if (!kinetic.data()) return;
                
                kinetic.run();
                kinetic.delete_kinetic_data_structure();

                std::cout << "finished building " << building_index << std::endl << std::endl;
            }

            void set_output(const JP_kinetic_propagation &kinetic, Building &building) const {

                building.polyhedrons.clear();
                for (JP_polyhedrons_iterator p_it = kinetic.partition->polyhedrons_begin(); p_it != kinetic.partition->polyhedrons_end(); ++p_it)
                    add_polyhedron(*p_it, building);

                add_polyhedron_facets(kinetic, building);
            }

            void add_polyhedron(const JP_polyhedron *input_polyhedron, Building &building) const {
                
                Polyhedrons &polyhedrons = building.polyhedrons;
                Polyhedron output_polyhedron;
                
                JP_sequences sequences_per_side; JP_conversions conversions;
                get_polyhedron_vertices(input_polyhedron, output_polyhedron.vertices, sequences_per_side, conversions);
                get_polyhedron_facets(sequences_per_side, conversions, output_polyhedron.facets);

                polyhedrons.push_back(output_polyhedron);
            }

            void get_polyhedron_vertices(const JP_polyhedron *polyhedron, Polyhedron_vertices &vertices, JP_sequences &sequences_per_side, JP_conversions &conversions) const {

                vertices.clear();
                sequences_per_side.clear();
                conversions.clear();

                JP_sequence_set vertices_used;
                for (JP_sides_iterator s_it = polyhedron->facets_begin(); s_it != polyhedron->facets_end(); ++s_it) {
                    const JP_facet *facet = s_it->first;

                    JP_sequence facet_vertices;
                    facet->get_circular_sequence_of_vertices(facet_vertices, !s_it->second);

                    for (JP_sequence_iterator v_it = facet_vertices.begin(); v_it != facet_vertices.end(); ++v_it) vertices_used.insert(*v_it);
                    sequences_per_side.push_back(facet_vertices);
                }

                for (JP_sequence_set_iterator v_it = vertices_used.begin(); v_it != vertices_used.end(); ++v_it) {
                    
                    const JP_vertex  *v = *v_it;
                    const JP_point_3 &p =  v->M;

                    const FT x = static_cast<FT>(CGAL::to_double(p.x()));
                    const FT y = static_cast<FT>(CGAL::to_double(p.y()));
                    const FT z = static_cast<FT>(CGAL::to_double(p.z()));
                    
                    const Polyhedron_vertex vertex = Polyhedron_vertex(x, y, z);

                    vertices.push_back(vertex);
                    conversions[v] = int(vertices.size()) - 1;
                }
            }

            void get_polyhedron_facets(const JP_sequences &sequences_per_side, const JP_conversions &conversions, Polyhedron_facets &facets) const {
                facets.clear();

                for (JP_sequences_iterator s_it = sequences_per_side.begin(); s_it != sequences_per_side.end(); ++s_it) {
                    const JP_sequence &sequence = *s_it;

                    Polyhedron_facet facet;
                    for (JP_sequence_iterator v_it = sequence.begin(); v_it != sequence.end(); ++v_it) facet.indices.push_back(conversions.at(*v_it));
                    facets.push_back(facet);
                }
            }

            void add_polyhedron_facets(const JP_kinetic_propagation &kinetic, Building &building) const {

                Polyhedrons &polyhedrons = building.polyhedron_facets;
                Polyhedron polyhedron;
                
                get_polyhedron_vertices_and_facets(kinetic, polyhedron.vertices, polyhedron.facets);
                polyhedrons.push_back(polyhedron);
            }

            void get_polyhedron_vertices_and_facets(const JP_kinetic_propagation &kinetic, Polyhedron_vertices &vertices, Polyhedron_facets &facets) const {

                JP_facet_vertices v;
                kinetic.partition->octree->get_all_vertices_sorted_by_identifier(v);
                
                CGAL_precondition(v.size() > 0);
                for (size_t i = 0; i < v.size(); ++i) {
                    const JP_point_3 &p = v[i]->M;

                    const FT x = static_cast<FT>(CGAL::to_double(p.x()));
                    const FT y = static_cast<FT>(CGAL::to_double(p.y()));
                    const FT z = static_cast<FT>(CGAL::to_double(p.z()));
                 
                    const Polyhedron_vertex vertex = Polyhedron_vertex(x, y, z);
                    vertices.push_back(vertex);
                }

                const auto &partiton_facets = kinetic.partition->facets;
                for (size_t i = 0; i < partiton_facets.size(); ++i) {
                    for (JP_facets_iterator f_it = partiton_facets[i].begin(); f_it != partiton_facets[i].end(); ++f_it) {
                        
                        const JP_facet* f = *f_it;
                        if (f->p < 6) continue;

                        JP_sequence facet_vertices;
                        f->get_circular_sequence_of_vertices(facet_vertices, true);

                        Polyhedron_facet facet;
                        for (JP_sequence_iterator v_it = facet_vertices.begin(); v_it != facet_vertices.end(); ++v_it) facet.indices.push_back((*v_it)->id);
                        facets.push_back(facet);
                    }
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_KINETIC_PARTITION_OUTPUT_CREATOR_STEP_8_H