#ifndef CGAL_LEVEL_OF_DETAIL_ROOF_PARTITIONING_TEST_H
#define CGAL_LEVEL_OF_DETAIL_ROOF_PARTITIONING_TEST_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_roof_partitioning_test {
            
        public:
            using Kernel    = InputKernel;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT        = typename Kernel::FT;
            using Point_2   = typename Kernel::Point_2;
            using Point_3   = typename Kernel::Point_3;

            using Wall  = typename Building::Wall;
            using Walls = typename Building::Walls;
            
            using Roof  = typename Building::Roof;
            using Roofs = typename Building::Roofs;

            using Buildings_iterator = typename Buildings::iterator;

            using Floor_faces = typename Building::Floor_faces;

            using Log = CGAL::LOD::Mylog;

            class Partitioning_data_structure {

            public:
                class Polygon_3 {

                public:
                    using Data = std::vector<Point_3>;

                    Data data;
                    bool propagate = false;
                };
                using Polygons_3 = std::vector<Polygon_3>;

                Polygons_3 polygons;
            };

            using Polygon_3  = typename Partitioning_data_structure::Polygon_3;
            using Polygons_3 = typename Partitioning_data_structure::Polygons_3;

            using Polygon_data = typename Polygon_3::Data;

            Level_of_detail_roof_partitioning_test(const FT ground_height, Buildings &buildings) :
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_big_value(FT(100000000000000)),
            m_z_scale(FT(10)),
            m_up_scale(FT(2)),
            m_down_scale(FT(1) / FT(2))
            { }

            void create() const {
                
                if (m_buildings.size() == 0)
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;
                    
					if (building.is_valid && building.shapes.size() != 0) 
                        process_building(building);
                } 
            }

        private:
            const FT m_ground_height;
            Buildings &m_buildings;
            
            const FT m_big_value;
            const FT m_z_scale;
            const FT m_up_scale;
            const FT m_down_scale;

            void process_building(Building &building) const {

                Partitioning_data_structure data;
                
                set_input(building, data);
                apply_partitioning(data);
            }

            void set_input(const Building &building, Partitioning_data_structure &data) const {

                data.polygons.clear();
                
                // set_ground(building, data);
                
                set_walls(building, data);
                set_roofs(building, data);

                Log logger;
                logger.export_partitioning_data_structure(data, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "7_data_structure_input", true);
            }

            void set_ground(const Building &building, Partitioning_data_structure &data) const {

                FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;

                const Floor_faces &faces = building.faces;
                for (size_t i = 0; i < faces.size(); ++i) {
                    
                    for (size_t j = 0; j < 3; ++j) {    
                        const Point_2 &p = faces[i]->vertex(j)->point();
                        
                        minx = CGAL::min(minx, p.x());
                        miny = CGAL::min(miny, p.y());

                        maxx = CGAL::max(maxx, p.x());
                        maxy = CGAL::max(maxy, p.y());
                    }
                }

                Polygon_data polygon_data(4);

                polygon_data[0] = Point_3(minx, miny, m_ground_height);
                polygon_data[1] = Point_3(maxx, miny, m_ground_height);
                polygon_data[2] = Point_3(maxx, maxy, m_ground_height);
                polygon_data[3] = Point_3(minx, maxy, m_ground_height);

                process_polygon(polygon_data, m_up_scale, FT(1), false, data.polygons);
            }

            void process_polygon(const Polygon_data &polygon_data, const FT scale, const FT z_extender, const bool propagate, Polygons_3 &polygons) const {

                if (polygon_data.size() == 0) 
                    return;
                
                Polygon_3 polygon;

                polygon.propagate = propagate;
                polygon.data      = polygon_data;

                scale_polygon(scale, z_extender, polygon);
                polygons.push_back(polygon);
            }

            void scale_polygon(const FT scale, const FT z_extender, Polygon_3 &polygon) const {

                Point_3 b;
                compute_barycentre(polygon, b);

                for (size_t i = 0; i < polygon.data.size(); ++i) {
                    Point_3 &p = polygon.data[i];

                    const FT x = (p.x() - b.x()) * scale + b.x();
                    const FT y = (p.y() - b.y()) * scale + b.y();
                    const FT z = (p.z() - b.z()) * scale * z_extender + b.z();

                    p = Point_3(x, y, z);
                }
            }

            void compute_barycentre(const Polygon_3 &polygon, Point_3 &b) const {

                CGAL_precondition(polygon.data.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < polygon.data.size(); ++i) {
                    const Point_3 &p = polygon.data[i];

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(polygon.data.size());
                y /= static_cast<FT>(polygon.data.size());
                z /= static_cast<FT>(polygon.data.size());

                b = Point_3(x, y, z);
            }

            void set_walls(const Building &building, Partitioning_data_structure &data) const {

                const Walls &walls = building.walls;
                for (size_t i = 0; i < walls.size(); ++i) {

                    const Wall &wall = walls[i];
                    process_polygon(wall.boundary, FT(1), m_z_scale, false, data.polygons);
                }
            }

            void set_roofs(const Building &building, Partitioning_data_structure &data) const {

                const Roofs &roofs = building.roofs;
                for (size_t i = 0; i < roofs.size(); ++i) {
                    
                    if (!roofs[i].is_valid)
                        continue;

                    const Roof &roof = roofs[i];
                    process_polygon(roof.boundary, m_up_scale, FT(1), true, data.polygons);
				}
            }

            void apply_partitioning(Partitioning_data_structure &data) const {



                Log logger;
                logger.export_partitioning_data_structure(data, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "7_data_structure_output", false);
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ROOF_PARTITIONING_TEST_H