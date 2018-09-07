#ifndef CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_BOX_STRATEGY_H
#define CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_BOX_STRATEGY_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Lod_2/Roofs_estimator/Level_of_detail_roofs_eigen_diagonalize_traits.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class ContainerInput, class BuildingInput>
		class Level_of_detail_roofs_estimator_box_strategy {
            
        public:
            typedef KernelTraits   Kernel;
            typedef ContainerInput Input;
            typedef BuildingInput  Building;

			typename Kernel::Compute_squared_length_3 squared_length_3;
			
			typename Kernel::Compute_scalar_product_2 dot_product_2;
			typename Kernel::Compute_scalar_product_3 dot_product_3;

			typename Kernel::Compute_determinant_2 			  cross_product_2;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using FT       = typename Kernel::FT;
            using Point_2  = typename Kernel::Point_2;
			using Point_3  = typename Kernel::Point_3;
            using Plane_3  = typename Kernel::Plane_3;
            using Line_2   = typename Kernel::Line_2;
            using Vector_2 = typename Kernel::Vector_2;
			using Vector_3 = typename Kernel::Vector_3;
            
			using Points = std::vector<Point_3>;

            using Local_kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

			using Local_ft      = typename Local_kernel::FT;
            using Local_point_2 = typename Local_kernel::Point_2;
            using Local_line_2  = typename Local_kernel::Line_2;

            using Diagonalize_traits = CGAL::LOD::Roofs_eigen_diagonalize_traits<Local_ft, 2>;

            using Roof = typename Building::Roof;

            Level_of_detail_roofs_estimator_box_strategy(const Input &input) : 
			m_input(input) 
			{ }

            void estimate_roof(Points &roof_points, const Plane_3 &plane, Building &building) const {
                if (roof_points.size() < 2) return;

				Vector_3 roof_normal;
				set_plane_normal(plane, roof_normal);

				Vector_3 ground_normal;
				set_ground_normal(ground_normal);

				FT angle_3d; Vector_3 axis;
				compute_angle_and_axis(roof_normal, ground_normal, angle_3d, axis);
				
				if (angle_3d != FT(0)) 
					rotate_points(angle_3d, axis, roof_points);

                Vector_2 roof_direction;
                estimate_roof_direction(roof_points, roof_direction);

                Vector_2 y_direction;
                set_y_direction(y_direction);

				FT angle_2d;
				compute_angle(roof_direction, y_direction, angle_2d);

				Point_2 barycentre;
                compute_barycentre(roof_points, barycentre);

				if (angle_2d != FT(0))
					rotate_points(angle_2d, barycentre, roof_points);

				Points boundary;
				compute_bounding_box(roof_points, boundary);

				rotate_points(-angle_2d, barycentre, boundary);
				rotate_points(-angle_3d, axis, boundary);
				
				if (!is_valid_boundary(boundary)) {
				
					boundary.clear();
					return;
				}

				create_roof(boundary, building);
            }

        private:
            const Input &m_input;

            void set_plane_normal(const Plane_3 &plane, Vector_3 &m) const {
				m = plane.orthogonal_vector();
			}

			void set_ground_normal(Vector_3 &n) const {

				const Plane_3 ground = Plane_3(FT(0), FT(0), FT(1), FT(0));				
				n = ground.orthogonal_vector();
			}

            void compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {
				
				const auto  cross = cross_product_3(m, n);
				const   FT length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const   FT    dot = dot_product_3(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
				if (angle == FT(0)) return;
				
				CGAL_precondition(length != FT(0));
				axis = cross / length;
			}

			void rotate_points(const FT angle, const Vector_3 &axis, Points &points) const {

				for (size_t i = 0; i < points.size(); ++i) {
					
					Point_3 &p = points[i];
					rotate_point(angle, axis, p);
				}
			}

			void rotate_point(const FT angle, const Vector_3 &axis, Point_3 &p) const {

				const double tmp_angle = CGAL::to_double(angle);

				const FT c = static_cast<FT>(std::cos(tmp_angle));
				const FT s = static_cast<FT>(std::sin(tmp_angle));

				const FT C = FT(1) - c;

				const FT x = axis.x();
				const FT y = axis.y();
				const FT z = axis.z();

				p = Point_3((x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
					  		(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
					  		(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
			}

			void project_points(Points &points) const {

				for (size_t i = 0; i < points.size(); ++i) {
					
					Point_3 &p = points[i];
					p = Point_3(p.x(), p.y(), FT(0));
				}
			}

            void estimate_roof_direction(const Points &points, Vector_2 &direction) const {

				Local_ft cx = Local_ft(0), cy = Local_ft(0);
                CGAL_precondition(points.size() > 1);

				std::vector<Local_point_2> tmp_points(points.size());
				for (size_t i = 0; i < points.size(); ++i) {
					
					const Point_3 &p = points[i];

					const Local_ft x = static_cast<Local_ft>(CGAL::to_double(p.x()));
					const Local_ft y = static_cast<Local_ft>(CGAL::to_double(p.y()));

					cx += x;
					cy += y;

					tmp_points[i] = Local_point_2(x, y);
				}

				cx /= static_cast<Local_ft>(points.size());
				cy /= static_cast<Local_ft>(points.size());

				Local_point_2 centroid = Local_point_2(cx, cy);

                Local_line_2 tmp_line;
				CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), tmp_line, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits());

				const Line_2 line = Line_2(
                Point_2(static_cast<FT>(tmp_line.point(0).x()), static_cast<FT>(tmp_line.point(0).y())), 
				Point_2(static_cast<FT>(tmp_line.point(1).x()), static_cast<FT>(tmp_line.point(1).y())));

                direction = line.to_vector();
            }

            void set_y_direction(Vector_2 &direction) const {
                direction = Vector_2(FT(0), FT(1));
            }

            void compute_angle(const Vector_2 &m, const Vector_2 &n, FT &angle) const {
				
				const FT cross = cross_product_2(m, n);
				const FT dot   = dot_product_2(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(cross), CGAL::to_double(dot)));
				if (angle == FT(0)) return;
			}

			void compute_barycentre(const Points &points, Point_2 &barycentre) const {

				FT bx = FT(0), by = FT(0);
				for (size_t i = 0; i < points.size(); ++i) {
							
					const Point_3 &p = points[i];

					bx += p.x();
					by += p.y();
				}

				bx /= static_cast<FT>(points.size());
                by /= static_cast<FT>(points.size());

                barycentre = Point_2(bx, by);
            }

			void rotate_points(const FT angle, const Point_2 &barycentre, Points &points) const {

				for (size_t i = 0; i < points.size(); ++i) {
					
					Point_3 &p = points[i];
					rotate_point(angle, barycentre, p);
				}
			}

            void rotate_point(const FT angle, const Point_2 &barycentre, Point_3 &p) const {

				FT x = p.x();
				FT y = p.y();

				x -= barycentre.x();
				y -= barycentre.y();

				p = Point_3(x, y, p.z());

                const double tmp_angle = CGAL::to_double(angle);

                const FT c = static_cast<FT>(std::cos(tmp_angle));
				const FT s = static_cast<FT>(std::sin(tmp_angle));

				x = p.x() * c - p.y() * s;
				y = p.y() * c + p.x() * s;

				x += barycentre.x();
				y += barycentre.y();

				p = Point_3(x, y, p.z());
			} 

			void compute_bounding_box(const Points &points, Points &bbox) const {

                const FT big_value = FT(100000000000000);

				FT minx =  big_value, miny =  big_value;
				FT maxx = -big_value, maxy = -big_value;

				FT z = FT(0);
				for (size_t i = 0; i < points.size(); ++i) {
					const Point_3 &p = points[i];

					minx = CGAL::min(minx, p.x());
					miny = CGAL::min(miny, p.y());

				    maxx = CGAL::max(maxx, p.x());
					maxy = CGAL::max(maxy, p.y());

					z += p.z();
				}
				z /= static_cast<FT>(points.size());

                bbox.clear();
				bbox.resize(4);

                bbox[0] = Point_3(minx, miny, z);
				bbox[1] = Point_3(maxx, miny, z);
				bbox[2] = Point_3(maxx, maxy, z);
				bbox[3] = Point_3(minx, maxy, z);
			}

			void create_roof(const Points &boundary, Building &building) const {

                Roof roof;

                roof.boundary = boundary;
				roof.is_valid = true;

                building.roofs.push_back(roof);
			}

			bool is_valid_boundary(const Points &boundary) const {

				if (std::isnan(CGAL::to_double(boundary[0].x())) ||
					std::isnan(CGAL::to_double(boundary[0].y())) ||
					std::isnan(CGAL::to_double(boundary[0].z()))  ) return false;

                if (std::isnan(CGAL::to_double(boundary[1].x())) ||
					std::isnan(CGAL::to_double(boundary[1].y())) ||
					std::isnan(CGAL::to_double(boundary[1].z()))  ) return false;

                if (std::isnan(CGAL::to_double(boundary[2].x())) ||
					std::isnan(CGAL::to_double(boundary[2].y())) ||
					std::isnan(CGAL::to_double(boundary[2].z()))  ) return false;

                if (std::isnan(CGAL::to_double(boundary[3].x())) ||
					std::isnan(CGAL::to_double(boundary[3].y())) ||
					std::isnan(CGAL::to_double(boundary[3].z()))  ) return false;

				if (boundary.size() < 3) return false;
				return true;
			}
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ROOFS_ESTIMATOR_BOX_STRATEGY_H