#ifndef CGAL_LEVEL_OF_DETAIL_TRIANGULATION_BASED_BOUNDARY_EXTRACTOR_H
#define CGAL_LEVEL_OF_DETAIL_TRIANGULATION_BASED_BOUNDARY_EXTRACTOR_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer>
		class Level_of_detail_triangulation_based_boundary_extractor {

		public:
            using Kernel = InputKernel;
            using Input  = InputContainer;

			using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

			using Projected_points = std::map<int, Point_2>;

			using Index   = int;
            using Indices = std::vector<Index>;

			Level_of_detail_triangulation_based_boundary_extractor(const Input &input, const Indices &boundary_indices, const Indices &interior_indices) :
            m_input(input),
            m_boundary_indices(boundary_indices),
            m_interior_indices(interior_indices)
            { }

            void extract(Projected_points &projected_points) const {

            }

		private:
			const Input   &m_input;
            const Indices &m_boundary_indices;
            const Indices &m_interior_indices;
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TRIANGULATION_BASED_BOUNDARY_EXTRACTOR_H