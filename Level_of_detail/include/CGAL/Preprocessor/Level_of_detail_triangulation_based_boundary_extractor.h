#ifndef CGAL_LEVEL_OF_DETAIL_TRIANGULATION_BASED_BOUNDARY_EXTRACTOR_H
#define CGAL_LEVEL_OF_DETAIL_TRIANGULATION_BASED_BOUNDARY_EXTRACTOR_H

// STL includes.
#include <map>
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Preprocessor/Level_of_detail_interior_boundary_extractor.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputCDT>
		class Level_of_detail_triangulation_based_boundary_extractor {

		public:
            using Kernel = InputKernel;
            using Input  = InputContainer;
            using CDT    = InputCDT;

			using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using Edge              = typename CDT::Edge;
            using Face_handle       = typename CDT::Face_handle;
            using Vertex_handle     = typename CDT::Vertex_handle;
            using Edges_iterator    = typename CDT::Finite_edges_iterator;
            using Vertices_iterator = typename CDT::Finite_vertices_iterator;
            using Vertex_circulator = typename CDT::Vertex_circulator;

			using Projected_points = std::map<int, Point_2>;

			using Index   = int;
            using Indices = std::vector<Index>;

            using Log = CGAL::LOD::Mylog;
            using Alpha_shapes_extractor = CGAL::LOD::Level_of_detail_interior_boundary_extractor<Kernel, Input, Indices>;

            using Heights = std::vector<FT>;

            typename Kernel::Compute_squared_distance_3 squared_distance_3;

            using Queue = std::queue<Vertex_handle>;

			Level_of_detail_triangulation_based_boundary_extractor(const Input &input, const Indices &boundary_indices, const Indices &interior_indices) :
            m_input(input),
            m_boundary_indices(boundary_indices),
            m_interior_indices(interior_indices),
            m_distance_threshold(FT(2)),
            m_use_std_method(false),
            m_check_uniqueness(true)
            { }

            void extract(Projected_points &projected_points) const {

                if (m_boundary_indices.size() != 0) return;

                CDT cdt;
                Queue queue;
                create_cdt_and_queue(projected_points, cdt, queue);

                if (m_use_std_method) extract_boundary_points_from_cdt_vertices(cdt, projected_points);
                else extract_boundary_points(cdt, queue, projected_points);
            }

		private:
			const Input   &m_input;
            const Indices &m_boundary_indices;
            const Indices &m_interior_indices;

            const FT m_distance_threshold;
            const bool m_use_std_method;

            const bool m_check_uniqueness;

            void create_cdt_and_queue(const Projected_points &projected_points, CDT &cdt, Queue &queue) const {
                
                cdt.clear();
                add_projected_points(projected_points, cdt, queue);

                add_boundary_points(cdt);
                add_interior_points(cdt);

                Log log;
                log.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "1_cdt");
            }

            void add_projected_points(const Projected_points &projected_points, CDT &cdt, Queue &queue) const {

                for (auto pit = projected_points.begin(); pit != projected_points.end(); ++pit) {
                    const auto &pair = *pit;

                    const Index index = pair.first;
                    const Point_2 &p  = pair.second;

                    queue.push(cdt.insert(p));
                    queue.front()->info().index = index;
                }
            }

            void add_boundary_points(CDT &cdt) const {

                for (size_t i = 0; i < m_boundary_indices.size(); ++i) {
                    
                    const Point_3 &p = m_input.point(m_boundary_indices[i]);
                    const Point_2  q = Point_2(p.x(), p.y());

                    Vertex_handle vh = cdt.insert(q);
                    
                    vh->info().height = p.z();
                    vh->info().index  = m_boundary_indices[i];
                }
            }

            void add_interior_points(CDT &cdt) const {

                for (size_t i = 0; i < m_interior_indices.size(); ++i) {
                    
                    const Point_3 &p = m_input.point(m_interior_indices[i]);
                    const Point_2  q = Point_2(p.x(), p.y());

                    Vertex_handle vh = cdt.insert(q);
                    
                    vh->info().height = p.z();
                    vh->info().index  = m_interior_indices[i];
                }
            }

            bool should_be_added(const Index index1, const Index index2) const {

                if (index1 < 0 || index1 >= m_input.number_of_points()) return false;
                if (index2 < 0 || index2 >= m_input.number_of_points()) return false;

                const Point_3 &p1 = m_input.point(index1);
                const Point_3 &p2 = m_input.point(index2);

                const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_3(p1, p2))));
                return distance > m_distance_threshold && distance < m_distance_threshold * FT(3);
            }

            void extract_boundary_points(const CDT &cdt, Queue &queue, Projected_points &projected_points) const {

                while (!queue.empty()) {
                    const Vertex_handle &curr = queue.front();
                    
                    Vertex_circulator vc = cdt.incident_vertices(curr);
                    add_boundary_point(curr, vc, queue, projected_points);
                }
            }

            void add_boundary_point(const Vertex_handle &curr, Vertex_circulator &vc, Queue &queue, Projected_points &projected_points) const {

                if (vc.is_empty()) {
                    
                    queue.pop();
                    return;
                }

                Vertex_circulator end = vc;
                do {
                    const Vertex_handle &vh = static_cast<Vertex_handle>(vc);
                    
                    if (should_be_added(curr->info().index, vh->info().index)) {
                        if (projected_points.find(vh->info().index) == projected_points.end()) {

                            projected_points[vh->info().index] = vh->point();
                            queue.push(vh);
                        }
                    }
                    ++vc;

                } while (vc != end);
                queue.pop();
            }

            void extract_boundary_points_from_cdt_vertices(const CDT &cdt, Projected_points &projected_points) const {

                for (Vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
                    const Vertex_handle &vh = static_cast<Vertex_handle>(vit);

                    Vertex_circulator vc = cdt.incident_vertices(vh);
                    add_boundary_point_from_cdt_vertices(vh, vc, projected_points);
                }
            }

            void add_boundary_point_from_cdt_vertices(const Vertex_handle &curr, Vertex_circulator &vc, Projected_points &projected_points) const {
                
                Indices indices;
                get_indices(vc, indices);

                if (!is_valid(curr, indices, projected_points)) return;

                for (size_t i = 0; i < indices.size(); ++i) {
                    if (should_be_added(curr->info().index, indices[i])) {
                     
                        projected_points[curr->info().index] = curr->point();
                        return;
                    }
                }   
            }

            void get_indices(Vertex_circulator &vc, Indices &indices) const {

                indices.clear();

                if (vc.is_empty()) return;
                Vertex_circulator end = vc;

                do {
                    const Vertex_handle &vh = static_cast<Vertex_handle>(vc);
                    indices.push_back(vh->info().index);
                    
                    ++vc;
                } while (vc != end);
            }

            bool is_valid(const Vertex_handle &curr, const Indices &indices, Projected_points &projected_points) const {

                if (projected_points.find(curr->info().index) != projected_points.end())
                    return false;

                if (indices.size() == 0) 
                    return false;

                if (m_check_uniqueness) {

                    for (size_t i = 0; i < indices.size(); ++i) {
                        const Index index = indices[i];

                        if (projected_points.find(index) != projected_points.end())
                            return false;
                    }
                }
                return true;
            }

            /*
            void extract_boundary_points_from_cdt_edges(const CDT &cdt, Projected_points &projected_points) const {

                projected_points.clear();
                for (Edges_iterator eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
                    
                    const Edge &edge = static_cast<Edge>(*eit);
                    add_boundary_point_from_cdt_edge(cdt, edge, projected_points);
                }
            }

            void add_boundary_point_from_cdt_edge(const CDT &cdt, const Edge &edge, Projected_points &projected_points) const {

                const Face_handle &fh1 = edge.first;
                const Index index1     = edge.second;

                const Face_handle &fh2 = fh1->neighbor(index1);
                const Index index2     = fh2->index(fh1);

                const Vertex_handle &vh1 = fh1->vertex(index1);
                const Vertex_handle &vh2 = fh2->vertex(index2);

                const FT height1 = vh1->info().height;
                const FT height2 = vh2->info().height;

                if (CGAL::abs(height1 - height2) > m_distance_threshold) {

                    const Index iv1 = (index1 + 1) % 3;
                    const Index iv2 = (index1 + 2) % 3;

                    const Vertex_handle &ev1 = fh1->vertex(iv1);
                    const Vertex_handle &ev2 = fh1->vertex(iv2);

                    projected_points[vh1->info().index] = vh1->point();
                    projected_points[vh2->info().index] = vh2->point();

                    projected_points[ev1->info().index] = ev1->point();
                    projected_points[ev2->info().index] = ev2->point();
                }
            }

            void apply_alpha_shapes(Projected_points &projected_points) const {

                Alpha_shapes_extractor extractor;
                extractor.set_alpha(FT(2));

                Indices stub, result;
                extractor.extract(m_input, stub, projected_points, result);

                projected_points.clear();
                for (size_t i= 0; i < result.size(); ++i) {
                    
                    const Point_3 &p = m_input.point(result[i]);
                    projected_points[result[i]] = Point_2(p.x(), p.y());
                }
            } */
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TRIANGULATION_BASED_BOUNDARY_EXTRACTOR_H