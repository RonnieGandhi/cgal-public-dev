#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_FACETS_CLEANER_3_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_FACETS_CLEANER_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>

// New CGAL includes.
#include <CGAL/Buildings/Utils/Level_of_detail_local_mesh_builder.h>
#include <CGAL/Region_growing/Level_of_detail_planar_region_growing.h>
#include <CGAL/Lod_2/Facet_regions_merger_3/Level_of_detail_facet_regions_merger_3.h>
#include <CGAL/Lod_2/Facets_based_region_growing_3/Level_of_detail_facets_based_region_growing_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_buildings_facets_cleaner_3 {
            
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

            using Mesh = CGAL::Polyhedron_3<Kernel>;
            
            using Mesh_facet_handle   = typename Mesh::Facet_const_handle;
            using Mesh_facets         = std::vector< std::vector<Mesh_facet_handle> >;
            using HDS                 = typename Mesh::HalfedgeDS;
            using Halfedge_handle     = typename Mesh::Halfedge_const_handle;
            using Facet_vertex_handle = typename Mesh::Facet::Vertex_const_handle;

            using Local_builder         = CGAL::LOD::Local_mesh_builder<Kernel, HDS, Building, Mesh_facet_handle>;
            using Planar_region_growing = CGAL::LOD::Level_of_detail_planar_region_growing<Kernel, Mesh, Mesh_facets>;

            using Building_region  = std::vector<Mesh_facet_handle>;
			using Building_regions = std::vector<Building_region>;
			using Regions 		   = std::vector<Building_regions>;

            using Clean_boundary = std::vector<Point_3>;
            using Clean_region   = std::vector<Clean_boundary>;
            using Clean_regions  = std::vector<Clean_region>;

            using Planar_region_merger        = CGAL::LOD::Level_of_detail_facet_regions_merger_3<Kernel, Building>;
            using Facets_based_region_growing = CGAL::LOD::Level_of_detail_facets_based_region_growing_3<Kernel>;

            Level_of_detail_buildings_facets_cleaner_3(Buildings &buildings) :
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(100000)),
            m_max_num_iters(50),
            m_use_global_conditions(true),
            m_use_triangulation_merging(true),
            m_use_all_facets(false),
            m_use_new_rg(false)
            { }

            void create_clean_facets() {
                
                // std::cout << m_use_all_facets << std::endl;

                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
            }

            void use_all_facets(const bool new_state) {
                m_use_all_facets = new_state;
            }

        private:
            Buildings &m_buildings;
            
            const FT     m_tolerance;
            const size_t m_max_num_iters;

            const bool m_use_global_conditions;
            const bool m_use_triangulation_merging;

            bool m_use_all_facets;
            const bool m_use_new_rg;
            
            void process_building(Building &building) const {
            
                Clean_facets &clean_facets = building.clean_facets;
                clean_facets.clear();
                
                building.is_clean = true;
                const Polyhedrons &polyhedrons = building.polyhedrons;
                
				for (size_t i = 0; i < polyhedrons.size(); ++i)
					process_polyhedron(i, polyhedrons, clean_facets);

                if (clean_facets.size() == 0) return;
                merge_clean_facets(clean_facets);
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

                    clean_facet.first.clear();
                    clean_facet.first.resize(facet.indices.size());

                    for (size_t j = 0; j < facet.indices.size(); ++j)
                        clean_facet.first[j] = vertices[facet.indices[j]];

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

                if (clean_facet.first.size() != facet.indices.size()) return false;

                size_t count = 0;
                for (size_t i = 0; i < clean_facet.first.size(); ++i) {
                    for (size_t j = 0; j < facet.indices.size(); ++j) {
                        
                        if (are_equal_points(clean_facet.first[i], vertices[facet.indices[j]])) {
                            
                            ++count;
                            break;
                        }
                    }
                }
                return count == clean_facet.first.size();
            }

            bool are_equal_points(const Point_3 &p, const Point_3 &q) const {
                return CGAL::abs(p.x() - q.x()) < m_tolerance && CGAL::abs(p.y() - q.y()) < m_tolerance && CGAL::abs(p.z() - q.z()) < m_tolerance;
            }

            void merge_clean_facets(Clean_facets &clean_facets) const {

                if (m_use_triangulation_merging) {
                    
                    merge_clean_facets_using_triangulation(clean_facets);
                    return;
                }

                Mesh mesh;

                Local_builder builder(clean_facets);
                mesh.delegate(builder);

                Mesh_facets facets;
                builder.get_facets(facets);

                Planar_region_growing planar_region_growing(facets);

                if (m_use_global_conditions) {
                    
                    planar_region_growing.use_global_conditions(true);
                    planar_region_growing.set_max_number_of_elements(1);
                }

				Regions regions;
				planar_region_growing.find_regions(regions);

                if (m_use_new_rg) {
                    
                    Clean_regions clean_regions;
                    Facets_based_region_growing region_growing(clean_facets);
                    region_growing.create_regions(clean_regions);
                    return;
                }

                set_merged_facets(regions[0], clean_facets);
            }

            void merge_clean_facets_using_triangulation(Clean_facets &clean_facets) const {

                Planar_region_merger planar_region_merger;
                planar_region_merger.use_all_facets(m_use_all_facets);

                planar_region_merger.merge(clean_facets);
            }

            void set_merged_facets(const Building_regions &regions, Clean_facets &clean_facets) const {
                
                clean_facets.clear();
                for (size_t i = 0; i < regions.size(); ++i)
                    create_clean_facet(regions[i], clean_facets);
            }

            void create_clean_facet(const Building_region &facets, Clean_facets &clean_facets) const {
                
                if (facets.size() == 1) {
                    
                    fall_back(facets, clean_facets);
                    return;
                }

                Clean_facet clean_facet;
                Halfedge_handle curr, end;
                size_t facet_index = 1000000;

                find_first_halfedge_handle(facets, curr, facet_index);
                CGAL_precondition(facet_index != 1000000);

                clean_facet.first.push_back(curr->vertex()->point());
                end = curr;

                curr = curr->next(); size_t iters = 0;
                while (curr != end) {

                    find_next_halfedge_handle(facets, facet_index, curr, clean_facet);

                    ++iters;
                    if (iters == m_max_num_iters) break;
                };
                
                if (iters == m_max_num_iters) {
                    
                    fall_back(facets, clean_facets);
                    return;
                }

                clean_facets.push_back(clean_facet);
            }

            void find_first_halfedge_handle(const Building_region &facets, Halfedge_handle &he, size_t &facet_index) const {

                CGAL_precondition(facets.size() != 0);

                for (size_t i = 0; i < facets.size(); ++i) {
                    facet_index = i;

                    const Mesh_facet_handle &fh = facets[facet_index];
                    
                    he = fh->halfedge();
                    if (is_unique(he, facet_index, facets)) return;

                    Halfedge_handle end = he;
                    he = he->next();

                    while (he != end) {

                        if (is_unique(he, facet_index, facets)) return;
                        he = he->next();
                    }
                }
            }

            bool is_unique(const Halfedge_handle &he, const size_t facet_index, const Building_region &facets) const {

                const Facet_vertex_handle &vh = he->vertex();
                const Point_3 &p = vh->point();

                Halfedge_handle stub;
                for (size_t i = 0; i < facets.size(); ++i) {
                    if (i == facet_index) continue;
                    
                    const Mesh_facet_handle &fh = facets[i];
                    if (contains(p, fh, stub)) return false;
                }
                return true;
            }

            bool contains(const Point_3 &p, const Mesh_facet_handle &fh, Halfedge_handle &next) const {

                Halfedge_handle curr = fh->halfedge();
                Halfedge_handle end  = curr;

                Facet_vertex_handle vh = curr->vertex();
                if (are_equal_points(p, vh->point())) {
                 
                    next = curr;
                    return true;
                }

                curr = curr->next();
                while (curr != end) {

                    vh = curr->vertex();
                    if (are_equal_points(p, vh->point())) {
                     
                        next = curr;
                        return true;
                    }

                    curr = curr->next();
                };
                return false;
            }

            void find_next_halfedge_handle(const Building_region &facets, size_t &facet_index, Halfedge_handle &curr, Clean_facet &clean_facet) const {

                const Point_3 &p = curr->vertex()->point();
                if (is_unique(curr, facet_index, facets)) {

                    clean_facet.first.push_back(p);
                    curr = curr->next();
                    
                    return;
                }

                Halfedge_handle next;
                for (size_t i = 0; i < facets.size(); ++i) {
                    if (i == facet_index) continue;
                    
                    const Mesh_facet_handle &fh = facets[i];
                    if (contains(p, fh, next)) {

                        clean_facet.first.push_back(p);
                        curr = next->next();

                        facet_index = i;
                        return;
                    }
                }
                return;
            }

            void fall_back(const Building_region &facets, Clean_facets &clean_facets) const {

                Clean_facet clean_facet;
                for (size_t i = 0; i < facets.size(); ++i) {

                    const Mesh_facet_handle &fh = facets[i];    
                    create_unique_clean_facet(fh, clean_facet);

                    clean_facets.push_back(clean_facet);
                }
            }

            void create_unique_clean_facet(const Mesh_facet_handle &fh, Clean_facet &clean_facet) const {

                clean_facet.first.clear();

                Halfedge_handle curr = fh->halfedge();
                Halfedge_handle end  = curr;

                Facet_vertex_handle vh = curr->vertex();
                clean_facet.first.push_back(vh->point());

                curr = curr->next();
                while (curr != end) {

                    vh = curr->vertex();
                    clean_facet.first.push_back(vh->point());

                    curr = curr->next();
                };
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_FACETS_CLEANER_3_H