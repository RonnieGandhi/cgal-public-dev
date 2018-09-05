#ifndef CGAL_LEVEL_OF_DETAIL_INSIDE_BUILDINGS_SELECTOR_STEP_0_H
#define CGAL_LEVEL_OF_DETAIL_INSIDE_BUILDINGS_SELECTOR_STEP_0_H

// STL includes.
#include <map>
#include <vector>
#include <unordered_set>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT, class InputBuildings>
		class Level_of_detail_inside_buildings_selector_step_0 {

        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputCDT       CDT;
            typedef InputBuildings Buildings;
			
            using Index   = int;
			using Indices = std::vector<Index>;

            using Vertex_handle = typename CDT::Vertex_handle;
            using Face_handle   = typename CDT::Face_handle;
            using Locate_type   = typename CDT::Locate_type;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;

            using Point_2    = typename Kernel::Point_2;
            using Triangle_2 = typename Kernel::Triangle_2;

            using Faces = std::vector<Face_handle>;

            using Building          = CGAL::LOD::Building<Kernel, CDT>;
            using Building_iterator = typename Buildings::iterator;

            Level_of_detail_inside_buildings_selector_step_0(const Input &input, const CDT &cdt, const Indices &indices) : 
            m_input(input), 
            m_cdt(cdt), 
            m_indices(indices),
            m_height_threshold(FT(1) / FT(1000000)) 
            { }

            void add_indices(Buildings &buildings) const {
                
                if (buildings.size() == 0) 
                    return;

                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    building.interior_indices.clear();
                    building.is_valid = true;
                }

                Locate_type locate_type;
				int locate_stub_index = -1;
                
                for (size_t i = 0; i < m_indices.size(); ++i) {
                    const Index &index = m_indices[i];

                    const Point_3 &point = m_input.point(index);
                    const Point_2 query  = Point_2(point.x(), point.y());

                    const Face_handle fh = m_cdt.locate(query, locate_type, locate_stub_index);
					if (locate_type == CDT::FACE || locate_type == CDT::EDGE || locate_type == CDT::VERTEX) {
                        
                        const int building_index = fh->info().bu;
                        if (building_index >= 0) buildings.at(building_index).interior_indices.push_back(index);
                    }
                }

                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    if (!is_valid_building(buildings, building)) building.is_valid = false; 
                    if (building.interior_indices.size() < 3)    building.is_valid = false;
                }
            }

        private:
            const Input   &m_input;
            const CDT     &m_cdt;
            const Indices &m_indices;

            const FT m_height_threshold;

			bool is_valid_building(const Buildings &buildings, Building &building) const {

				const FT height = building.height;
				if (height < m_height_threshold) return false;

				const auto &faces = building.faces;
				if (faces.size() < 2) {
				
					for (std::unordered_set<int>::const_iterator nit = building.neighbours.begin(); nit != building.neighbours.end(); ++nit) {
						if (is_valid_local_building(buildings.at(*nit))) {

							building.height = buildings.at(*nit).height;
							building.color  = buildings.at(*nit).color;

							return true;
						}
					}
					return false;
				}
				return true;
			}

			bool is_valid_local_building(const Building &building) const {

				const FT height = building.height;
				if (height < m_height_threshold) return false;

				const auto &faces = building.faces;
				if (faces.size() < 2) return false;

				return true;
			}

            void add_indices_to_building(Building &building) const {
                const Faces &building_faces = building.faces;
                
                for (size_t i = 0; i < building_faces.size(); ++i)
                    add_face_indices(building_faces[i], building);
            }

            void add_face_indices(const Face_handle &fh, Building &building) const {

                for (size_t i = 0; i < m_indices.size(); ++i) {
                    const Index index = m_indices[i];

                    const Point_3 &query = m_input.point(index);
                    if (belongs_to_face(query, fh)) building.interior_indices.push_back(index); 
                }
            }

            bool belongs_to_face(const Point_3 &query, const Face_handle &fh) const {

				const Point_2 &p1 = fh->vertex(0)->point();
				const Point_2 &p2 = fh->vertex(1)->point();
				const Point_2 &p3 = fh->vertex(2)->point();

				const Triangle_2 triangle = Triangle_2(p1, p2, p3);
				const Point_2 new_query   = Point_2(query.x(), query.y());

				if (triangle.has_on_bounded_side(new_query) || triangle.has_on_boundary(new_query)) return true;
				return false;
			}
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_INSIDE_BUILDINGS_SELECTOR_STEP_0_H