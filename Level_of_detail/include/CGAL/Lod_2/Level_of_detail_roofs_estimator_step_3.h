#ifndef CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_STEP_3_H
#define CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_STEP_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Lod_2/Roofs_estimator/Level_of_detail_roofs_diagonalize_traits.h>
#include <CGAL/Lod_2/Roofs_estimator/Level_of_detail_roofs_estimator_box_strategy.h>

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

            using Local_kernel       = CGAL::Simple_cartesian<double>;
            using Diagonalize_traits = CGAL::LOD::Eigen_roofs_diagonalize_traits<double, 3>;

			using Point_3ft = typename Local_kernel::Point_3;
            using Plane_3ft = typename Local_kernel::Plane_3;

            using Points_3ft = std::vector<Point_3ft>;

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
                building.clear_planes();

                const Shapes &shapes = building.shapes;
                if (shapes.size() == 0) {
                 
                    building.is_valid = false;
                    return;
                }

                building.clear_roofs();
				for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    process_roof(indices, building);
                }
            }

            void process_roof(const Indices &indices, Building &building) const {
                
                CGAL_precondition(indices.size() > 2);

                Plane_3 plane;
                fit_plane_to_roof_points(indices, plane);
                building.planes.push_back(plane);

                Points_3 points;   
                project_points_onto_plane(indices, plane, points);
                m_strategy.estimate_roof(points, plane, building);
            }

            void fit_plane_to_roof_points(const Indices &indices, Plane_3 &plane) const {
                
                Point_3ft centroid;
                Points_3ft points;
                set_points_and_centroid(indices, points, centroid);

                Plane_3ft tmp_plane;
				CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), tmp_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits());
				plane = Plane_3(static_cast<FT>(tmp_plane.a()), static_cast<FT>(tmp_plane.b()), static_cast<FT>(tmp_plane.c()), static_cast<FT>(tmp_plane.d()));
            }

            void set_points_and_centroid(const Indices &indices, Points_3ft &points, Point_3ft &centroid) const {
                CGAL_precondition(indices.size() > 2);

                points.clear();
                points.resize(indices.size());

                double bx = 0.0, by = 0.0, bz = 0.0;
				for (size_t i = 0; i < indices.size(); ++i) {

					const Point_3 &p = m_input.point(indices[i]);

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());
					const double z = CGAL::to_double(p.z());

					points[i] = Point_3ft(x, y, z);

                    bx += x;
                    by += y;
                    bz += z;
				}

                bx /= static_cast<double>(indices.size());
                by /= static_cast<double>(indices.size());
                bz /= static_cast<double>(indices.size());

                centroid = Point_3ft(bx, by, bz);
            }

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
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_STEP_3_H