#ifndef CGAL_LEVEL_OF_DETAIL_LOD0_LOD1_H
#define CGAL_LEVEL_OF_DETAIL_LOD0_LOD1_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Lod_2/Facet_regions_merger_3/Level_of_detail_facet_regions_merger_3.h>
#include <CGAL/Lod_2/Facets_based_region_growing_3/Level_of_detail_facets_based_region_growing_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputGround, class InputCDT, class InputBuilding, class InputBuildings>
		class Level_of_detail_lod0_lod1 {
            
        public:
            using Kernel    = InputKernel;
            using Ground    = InputGround;
            using CDT       = InputCDT;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using Buildings_iterator = typename Buildings::iterator;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_3 = typename Kernel::Triangle_3;

            using Log = CGAL::LOD::Mylog;

            using Facet_regions_merger        = Level_of_detail_facet_regions_merger_3<Kernel, Building>;
            using Facets_based_region_growing = CGAL::LOD::Level_of_detail_facets_based_region_growing_3<Kernel>;

            using Input_facet    = typename Facets_based_region_growing::Input_facet;
            using Input_facets   = typename Facets_based_region_growing::Input_facets;

            using Output_regions = typename Facets_based_region_growing::Output_regions;

            using Boundary = typename Building::Boundary;

            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;

            using Roof  = typename Building::Roof;
            using Roofs = typename Building::Roofs;

            using Face_handle = typename Building::Face_handle;
            using Floor_faces = typename Building::Floor_faces;

            using Wall_boundary = typename Wall::Wall_boundary;
            using Roof_boundary = typename Roof::Roof_boundary;

            using Output_facet  = typename Building::Region_facet;
            using Output_facets = std::vector<Output_facet>;

            using Color = CGAL::Color;
			
            Level_of_detail_lod0_lod1(const Ground &ground, const FT ground_height, const CDT &cdt, Buildings &buildings) :
            m_ground(ground),
            m_tolerance(FT(1) / FT(1000)),
            m_ground_height(ground_height),
            m_cdt(cdt),
            m_buildings(buildings)
            { }

            void create() {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
                save_buildings();
            }

        private:
            const CDT    &m_cdt;
            const Ground &m_ground;
            
            const FT m_tolerance;
            const FT m_ground_height;
            
            Buildings &m_buildings;

            void process_building(Building &building) const {

                if (building.faces.size() == 1) {
                 
                    building.is_valid = false;
                    return;
                }

                create_walls(building);
                create_roofs(building);
            }

            void create_walls(Building &building) const {
                
                detect_walls(building);
                merge_walls(building);
            }

            void detect_walls(Building &building) const {

                const Boundary &boundary       = building.boundaries[0];
                Output_regions &output_regions = building.output_regions;
                
                Input_facets input_facets;
                create_wall_input_facets(boundary, m_ground_height, m_ground_height + building.height, input_facets);

                Facets_based_region_growing region_growing(input_facets);
                region_growing.create_regions(output_regions);
            }

            void create_wall_input_facets(const Boundary &boundary, const FT height_floor, const FT height_roof, Input_facets &input_facets) const {
                
                input_facets.clear();
				for (size_t i = 0; i < boundary.size(); i += 2) {

					const size_t ip = i + 1;
					add_wall_input_facet(boundary[i]->point(), boundary[ip]->point(), height_floor, height_roof, input_facets);
				}
            }

			void add_wall_input_facet(const Point_2 &a, const Point_2 &b, const FT height_floor, const FT height_roof, Input_facets &input_facets) const {

                Input_facet input_facet(4);

				input_facet[0] = Point_3(a.x(), a.y(), height_floor);
				input_facet[1] = Point_3(b.x(), b.y(), height_floor);

				input_facet[2] = Point_3(b.x(), b.y(), height_roof);
				input_facet[3] = Point_3(a.x(), a.y(), height_roof);

                input_facets.push_back(input_facet);
			}

            void merge_walls(Building &building) const {

                Walls &walls = building.walls;
                Output_facets output_facets;

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(output_facets);

                set_final_wall_facets(output_facets, walls);
            }

            void set_final_wall_facets(const Output_facets &output_facets, Walls &walls) const {

                walls.clear();
                walls.resize(output_facets.size());

                for (size_t i = 0; i < output_facets.size(); ++i)
                    set_final_wall_facet(output_facets[i], walls[i]);
            }

            void set_final_wall_facet(const Output_facet &output_facet, Wall &wall) const {

                Wall_boundary &boundary = wall.boundary;
                boundary.clear();

                const size_t n = output_facet.size();
                for (size_t i = 0; i < n; ++i) {
                    
                    const size_t im = (i + n - 1) % n;
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

            void create_roofs(Building &building) const {
                
                detect_roofs(building);
                merge_roofs(building);
            }

            void detect_roofs(Building &building)const {

                const Floor_faces &faces = building.faces;
                Output_regions &output_regions = building.output_regions;

                Input_facets input_facets;
                create_roof_input_facets(faces, m_ground_height + building.height, input_facets);
                create_roof_output_regions(input_facets, output_regions);

                // Facets_based_region_growing region_growing(input_facets);
                // region_growing.create_regions(output_regions);
            }

            void create_roof_input_facets(const Floor_faces &faces, const FT height, Input_facets &input_facets) const {
                
                input_facets.clear();
				for (size_t i = 0; i < faces.size(); ++i)
					add_roof_input_facet(faces[i], height, input_facets);
            }

            void add_roof_input_facet(const Face_handle &fh, const FT height, Input_facets &input_facets) const {

                const Point_2 &a = fh->vertex(0)->point();
                const Point_2 &b = fh->vertex(1)->point();
                const Point_2 &c = fh->vertex(2)->point();

                Input_facet input_facet(3);

                input_facet[0] = Point_3(a.x(), a.y(), height);
                input_facet[1] = Point_3(b.x(), b.y(), height);
                input_facet[2] = Point_3(c.x(), c.y(), height);

                input_facets.push_back(input_facet);
            }

            void create_roof_output_regions(const Input_facets &input_facets, Output_regions &output_regions) const {

                output_regions.clear();
                output_regions.push_back(input_facets);
            }

            void merge_roofs(Building &building) const {

                Roofs &roofs = building.roofs;
                Output_facets output_facets;

                Facet_regions_merger facet_regions_merger(building.output_regions);
                facet_regions_merger.merge(output_facets);

                set_final_roof_facets(output_facets, roofs);
            }

            void set_final_roof_facets(const Output_facets &output_facets, Roofs &roofs) const {

                roofs.clear();
                roofs.resize(output_facets.size());

                for (size_t i = 0; i < output_facets.size(); ++i)
                    set_final_roof_facet(output_facets[i], roofs[i]);
            }

            void set_final_roof_facet(const Output_facet &output_facet, Roof &roof) const {

                Roof_boundary &boundary = roof.boundary;
                boundary.clear();

                const size_t n = output_facet.size();
                for (size_t i = 0; i < n; ++i) {
                    
                    const size_t im = (i + n - 1) % n;
                    const size_t ip = (i + 1) % n;

                    const Point_3 &p1 = output_facet[im];
                    const Point_3 &p2 = output_facet[i];
                    const Point_3 &p3 = output_facet[ip];

                    if (!are_collinear(p1, p2, p3)) boundary.push_back(p2);
                }
            }

            void save_buildings() const {

                Log exporter; 
                exporter.save_building_walls(m_buildings, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "11_walls", true);
                exporter.save_building_roofs_without_faces(m_buildings, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "12_roofs", true);

                const Color ground_color = Color(186, 189, 182);
                const Color wall_color   = Color(255, 255, 255);
                const Color roof_color   = Color(245, 121, 0);

                Ground ground;
                fix_ground(ground);

                exporter.save_lod0(m_buildings, ground, m_ground_height, ground_color, roof_color, "LOD0");
                exporter.save_lod1(m_buildings, ground, ground_color, wall_color, roof_color, "LOD1");
            }

            void fix_ground(Ground &ground) const {

                ground.clear();
                ground.resize(4);

                for (size_t i = 0; i < m_ground.size(); ++i) {
                    
                    const Point_3 &p = m_ground[i];
                    ground[i] = Point_3(p.x(), p.y(), m_ground_height);
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LOD0_LOD1_H