#ifndef CGAL_LEVEL_OF_DETAIL_REGION_GROWING_STEP_1_H
#define CGAL_LEVEL_OF_DETAIL_REGION_GROWING_STEP_1_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <unordered_set>

// CGAL includes.
#include <CGAL/Point_set_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Lod_2/Region_growing_simon_3/Plane.h>
#include <CGAL/Lod_2/Region_growing_simon_3/Shape_base.h>
#include <CGAL/Lod_2/Region_growing_simon_3/Region_growing.h>
#include <CGAL/Lod_2/Region_growing_simon_3/Shape_detection_traits.h>
#include <CGAL/Lod_2/Region_growing_simon_3/regularize_planes.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT, class InputBuildings>
		class Level_of_detail_region_growing_step_1 {

		public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputCDT       CDT;
            typedef InputBuildings Buildings;

            using Index   = int;
			using Indices = std::vector<Index>;

            using FT       = typename Kernel::FT;
            using Point_3  = typename Kernel::Point_3;
            using Vector_3 = typename Kernel::Vector_3;
            using Plane_3  = typename Kernel::Plane_3;

            using Vertex_handle   = typename CDT::Vertex_handle;
            using Face_handle     = typename CDT::Face_handle;

            using Building          = CGAL::LOD::Building<Kernel, CDT>;
            using Building_iterator = typename Buildings::iterator;

            using Log = CGAL::LOD::Mylog;

            using Final_planes = typename Building::Planes;

            using Local_kernel = CGAL::Simple_cartesian<double>;
			using Point_3ft    = Local_kernel::Point_3;
			using Vector_3ft   = Local_kernel::Vector_3;

            using Point_set  = CGAL::Point_set_3<Point_3ft>;
            using Point_map  = typename Point_set::Point_map;
            using Vector_map = typename Point_set::Vector_map;

            using Traits             = CGAL::Shape_detection_simon_3::Shape_detection_traits<Local_kernel, Point_set, Point_map, Vector_map>;
            using Region_growing     = CGAL::Shape_detection_simon_3::Region_growing<Traits>;
            using Plane              = CGAL::Shape_detection_simon_3::Plane<Traits>;
            using Plane_map          = CGAL::Shape_detection_simon_3::Plane_map<Traits>;
            using Point_to_plane_map = CGAL::Shape_detection_simon_3::Point_to_shape_index_map<Traits>;

            using Parameters     = typename CGAL::Shape_detection_simon_3::Region_growing<Traits>::Parameters;
            
            using Shape_iterator = typename Region_growing::Shape_range::iterator;
            using Shape          = typename Region_growing::Shape;
            using Planes         = typename Region_growing::Plane_range;
            
            using Index_iterator = std::vector<size_t>::const_iterator;

            Level_of_detail_region_growing_step_1(const Input &input) :
            m_input(input), 
            m_epsilon(-FT(1)),
			m_cluster_epsilon(-FT(1)),
			m_normal_threshold(-FT(1)),
			m_min_points(0),
            m_regularize(false),
            m_regularize_parallelism(true),
            m_regularize_orthogonality(true),
            m_regularize_coplanarity(false),
            m_regularize_axis_symmetry(true),
            m_regularization_angle(FT(25))
            { }

            void detect(Buildings &buildings) const {
                
                if (buildings.size() == 0) return;
                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    
                    Building &building = (*bit).second;
                    grow_regions(building);
                }
            }

            void regularize(const bool new_state) {
                m_regularize = new_state;
            }

            void set_regularization_angle(const FT new_value) {
                m_regularization_angle = new_value;
            }

            void set_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_epsilon = new_value;
			}

			void set_cluster_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_cluster_epsilon = new_value;
			}

			void set_normal_threshold(const FT new_value) {

				assert(new_value > FT(0) && new_value < FT(1));
				m_normal_threshold = new_value;
			}

			void set_minimum_shape_points(const size_t new_value) {

				assert(new_value > 0);
				m_min_points = new_value;
			}

        private:
            const Input &m_input;

            FT m_epsilon;
            FT m_cluster_epsilon;
            FT m_normal_threshold;
            size_t m_min_points;

            bool m_regularize;

            const bool m_regularize_parallelism;
            const bool m_regularize_orthogonality;
            const bool m_regularize_coplanarity;
            const bool m_regularize_axis_symmetry;

            FT m_regularization_angle;

            void grow_regions(Building &building) const {

                const Indices &indices = building.interior_indices;
                if (indices.size() == 0) return;

                Point_set points;
                set_points(indices, points);

                apply_3d_region_growing(indices, points, building);
            }

            void set_points(const Indices &indices, Point_set &points) const {
                
                points.clear();
                points.add_normal_map();

                for (size_t i = 0; i < indices.size(); ++i) {
                    const Index index = indices[i];

                    const Point_3  &point  = m_input.point(index);
                    const Vector_3 &normal = m_input.normal(index);

                    double x = CGAL::to_double(point.x());
                    double y = CGAL::to_double(point.y());
                    double z = CGAL::to_double(point.z());

                    const Point_3ft tmp_point = Point_3ft(x, y, z);

                    x = CGAL::to_double(normal.x());
                    y = CGAL::to_double(normal.y());
                    z = CGAL::to_double(normal.z());

                    const Vector_3ft tmp_normal = Vector_3ft(x, y, z);

                    points.insert(tmp_point, tmp_normal);
                }
            }

            void apply_3d_region_growing(const Indices &indices, Point_set &points, Building &building) const {

                Region_growing region_growing;

                region_growing.set_input(points, points.point_map(), points.normal_map());
                region_growing. template add_shape_factory<Plane>();

                assert(m_epsilon          > FT(0));
                assert(m_cluster_epsilon  > FT(0));
                assert(m_normal_threshold > FT(0) && m_normal_threshold < FT(1));
                assert(m_min_points       > 0);

                Parameters parameters;
                set_parameters(parameters);

                region_growing.detect(parameters);
                
                Final_planes &final_planes = building.planes;
                final_planes.clear();

                const Planes &planes = region_growing.planes();
                if (m_regularize) regularize_shapes(points, planes);
                
                set_final_planes(planes, final_planes);
                set_shapes_to_building(region_growing, indices, building);
            }

            void set_parameters(Parameters &parameters) const {

                parameters.epsilon          = CGAL::to_double(m_epsilon / FT(4));
                parameters.cluster_epsilon  = CGAL::to_double(m_cluster_epsilon);
                parameters.normal_threshold = CGAL::to_double(m_normal_threshold);
                parameters.min_points       = m_min_points * 6;
            }

            void regularize_shapes(const Point_set &points, const Planes &planes) const {

                CGAL::Shape_detection_simon_3::regularize_planes(points, points.point_map(), planes, Plane_map(), Point_to_plane_map(points, planes), 
                m_regularize_parallelism, m_regularize_orthogonality, m_regularize_coplanarity, m_regularize_axis_symmetry, CGAL::to_double(m_regularization_angle));
            }

            void set_final_planes(const Planes &planes, Final_planes &final_planes) const {

                final_planes.clear();
                final_planes.resize(planes.size());

                for (size_t i = 0; i < planes.size(); ++i) {
                    
                    const auto &plane = get(Plane_map(), *(planes.begin() + i));
                    final_planes[i] = Plane_3(static_cast<FT>(plane.a()), static_cast<FT>(plane.b()), static_cast<FT>(plane.c()), static_cast<FT>(plane.d()));
                }
            }

            void set_shapes_to_building(const Region_growing &region_growing, const Indices &indices, Building &building) const {
                
                building.clear_shapes();
                const size_t number_of_shapes = static_cast<size_t>(region_growing.shapes().end() - region_growing.shapes().begin());

                building.shapes.resize(number_of_shapes);
                Shape_iterator sit = region_growing.shapes().begin();

                size_t count = 0;
                while (sit != region_growing.shapes().end()) {
                    boost::shared_ptr<Shape> shape = *sit;

                    const size_t number_of_indices = static_cast<size_t>(shape->indices_of_assigned_points().end() - shape->indices_of_assigned_points().begin());
                    building.shapes[count].resize(number_of_indices);

                    Index_iterator index_it = shape->indices_of_assigned_points().begin(); size_t i = 0;
                    while (index_it != shape->indices_of_assigned_points().end()) {
                        const auto index = *index_it;

                        assert(index >= 0);
                        assert(index < indices.size());

                        building.shapes[count][i] = indices[index];
                           
                        ++index_it;
                        ++i;
                    }
                    ++sit;
                    ++count;
                }
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_REGION_GROWING_STEP_1_H