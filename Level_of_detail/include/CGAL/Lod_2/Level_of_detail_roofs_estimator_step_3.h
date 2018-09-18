#ifndef CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_STEP_3_H
#define CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_STEP_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// New CGAL includes.
#include <CGAL/Lod_2/Roofs_estimator/Level_of_detail_roofs_estimator_box_strategy.h>
#include <CGAL/Lod_2/Roofs_estimator/Level_of_detail_roofs_eigen_diagonalize_traits.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_roofs_estimator_step_3 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputContainer Input;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            using Buildings_iterator = typename Buildings::iterator;

            using Index   = int;
			using Indices = std::vector<Index>;
            using Shapes  = typename Building::Shapes;
            
            using Points_3 = std::vector<Point_3>;

            using Local_kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

			using Local_ft      = typename Local_kernel::FT;
            using Local_point_3 = typename Local_kernel::Point_3;
            using Local_plane_3 = typename Local_kernel::Plane_3;

            using Diagonalize_traits = CGAL::LOD::Roofs_eigen_diagonalize_traits<Local_ft, 3>;

            using Local_points_3 = std::vector<Local_point_3>;

            using Roof_estimation_strategy = CGAL::LOD::Level_of_detail_roofs_estimator_box_strategy<Kernel, Input, Building>;

            Level_of_detail_roofs_estimator_step_3(const Input &input, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_strategy(input)
            { }

            void estimate() {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;
            
            Roof_estimation_strategy m_strategy;

            void process_building(Building &building) const {

                if (building.planes.size() != building.shapes.size())
                    building.planes.clear();

                const Shapes &shapes = building.shapes;
                if (shapes.size() == 0) {
                 
                    building.is_valid = false;
                    return;
                }

                building.clear_roofs();
				for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    process_roof(indices, i, building);
                }
            }

            void process_roof(const Indices &indices, const size_t plane_index, Building &building) const {
                
                CGAL_precondition(indices.size() > 2);
                Plane_3 plane;

                // if (building.planes.size() == building.shapes.size()) plane = building.planes[plane_index];
                // else {
                    
                //     fit_plane_to_roof_points(indices, plane);
                //     building.planes.push_back(plane);
                // }

                plane = building.planes[plane_index];

                Points_3 points;   
                project_points_onto_plane(indices, plane, points);
                m_strategy.estimate_roof(points, plane, building);
            }

            /*
            void fit_plane_to_roof_points(const Indices &indices, Plane_3 &plane) const {
                
                Local_point_3  centroid;
                Local_points_3 points;
                set_points_and_centroid(indices, points, centroid);

                Local_plane_3 tmp_plane;
				CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), tmp_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits());
				plane = Plane_3(static_cast<FT>(tmp_plane.a()), static_cast<FT>(tmp_plane.b()), static_cast<FT>(tmp_plane.c()), static_cast<FT>(tmp_plane.d()));
            }

            void set_points_and_centroid(const Indices &indices, Local_points_3 &points, Local_point_3 &centroid) const {
                CGAL_precondition(indices.size() > 2);

                points.clear();
                points.resize(indices.size());

                Local_ft cx = Local_ft(0), cy = Local_ft(0), cz = Local_ft(0);
				for (size_t i = 0; i < indices.size(); ++i) {

					const Point_3 &p = m_input.point(indices[i]);

					const Local_ft x = static_cast<Local_ft>(CGAL::to_double(p.x()));
					const Local_ft y = static_cast<Local_ft>(CGAL::to_double(p.y()));
					const Local_ft z = static_cast<Local_ft>(CGAL::to_double(p.z()));

					points[i] = Local_point_3(x, y, z);

                    cx += x;
                    cy += y;
                    cz += z;
				}

                cx /= static_cast<Local_ft>(indices.size());
                cy /= static_cast<Local_ft>(indices.size());
                cz /= static_cast<Local_ft>(indices.size());

                centroid = Local_point_3(cx, cy, cz);
            } */

            void project_points_onto_plane(const Indices &indices, const Plane_3 &plane, Points_3 &points) const {
                CGAL_precondition(indices.size() > 2);

                points.clear();
                points.resize(indices.size());

                for (size_t i = 0; i < indices.size(); ++i) {			
					
                    const Point_3 &p = m_input.point(indices[i]);
					points[i] = plane.projection(p);
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_STEP_3_H