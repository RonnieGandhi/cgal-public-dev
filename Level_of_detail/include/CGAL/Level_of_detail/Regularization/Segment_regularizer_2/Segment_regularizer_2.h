#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Segment_regularizer_parameters.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Angles_regularizer.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Ordinates_regularizer.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment_property_map.h>

#include "../../../../../test/Level_of_detail/include/debugging/Mydebugger.h"

namespace CGAL {

	namespace Level_of_detail {

        template<class KernelTraits>
        class Level_of_detail_segment_regularizer_2 {

        public:
            typedef KernelTraits Kernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Line    = typename Kernel::Line_2;

            using Regular_segment     = CGAL::Level_of_detail::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments    = std::vector<Regular_segment *>;
            using In_regular_segments = std::vector<Regular_segment>;

            using RegularMap   = CGAL::Level_of_detail::Level_of_detail_segment_regularizer_regular_segment_property_map<Regular_segment, Segment>;
            using RegularRange = Regular_segments;
            
            using Angle_regularizer    = CGAL::Level_of_detail::Level_of_detail_segment_regularizer_for_angles<Kernel>;
            using Ordinate_regularizer = CGAL::Level_of_detail::Level_of_detail_segment_regularizer_for_ordinates<Kernel>;

            using Logger     = CGAL::Level_of_detail::Mydebugger;
            using Parameters = CGAL::Level_of_detail::Level_of_detail_segment_regularizer_parameters<Kernel>;

            using Tree = typename Angle_regularizer::Tree;

            Level_of_detail_segment_regularizer_2() : 
            m_silent(false), m_optimize_angles(true), m_optimize_ordinates(true) { }

            ~Level_of_detail_segment_regularizer_2() {
                m_input_segments.clear();
            }

            template<typename SegmentRange, typename SegmentMap>
            void regularize(SegmentRange &input_segments, SegmentMap segment_map) {
                if (input_segments.size() == 0) return;

                // Copy input segments in the internal storage.
                copy_input_segments(input_segments, segment_map);
                
                // Make segments parallel and orthogonal.
                if (m_optimize_angles) 
                    apply_angle_regularization();

                // Make segments collinear.
                if (m_optimize_angles && m_optimize_ordinates) 
                    apply_ordinate_regularization();

                // Update orientations of input segments.
                update_input_segments(input_segments, segment_map);
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

            void optimize_angles(const bool new_state) {
                m_optimize_angles = new_state;
            }

            void optimize_ordinates(const bool new_state) {
                m_optimize_ordinates = new_state;
            }

            // LOD function, can be removed.
            void get_lines_from_segments(const std::vector<Segment> &segments, std::vector<Line> &lines) const {
                
                lines.clear();
                lines.resize(segments.size());

                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point &source = segments[i].source();
                    const Point &target = segments[i].target();

                    lines[i] = Line(source, target);
                }
            }

            void set_max_angle_in_degrees(const FT new_value) {
                assert(new_value > FT(0));
                m_parameters.set_max_angle_in_degrees(new_value);
            }

            void set_max_difference_in_meters(const FT new_value) {
                assert(new_value > FT(0));
                m_parameters.set_max_difference_in_meters(new_value);
            }

        private:
            bool   m_silent;
            Logger m_logger;
            
            Regular_segments    m_input_segments;
            In_regular_segments m_in_segments;
            Parameters          m_parameters;

            bool m_optimize_angles;
            bool m_optimize_ordinates;

            std::shared_ptr<Angle_regularizer>    m_angle_regularizer_ptr;
            std::shared_ptr<Ordinate_regularizer> m_ordinate_regularizer_ptr;
            
            template<typename SegmentRange, typename SegmentMap>
            void copy_input_segments(const SegmentRange &input_segments, const SegmentMap &segment_map) {
                
                using Segment_iterator = typename SegmentRange::const_iterator;
                assert(input_segments.size() > 0);

                m_in_segments.clear();
                m_in_segments.resize(input_segments.size());

                size_t i = 0;
                for (Segment_iterator sit = input_segments.begin(); sit != input_segments.end(); ++sit, ++i) {
                    
                    const Segment &segment = get(segment_map, *sit);
                    m_in_segments[i] = Regular_segment(i, segment);
                }

                assert(i == m_in_segments.size());

                m_input_segments.clear();
                m_input_segments.resize(m_in_segments.size());

                for (size_t i = 0; i < m_in_segments.size(); ++i)
                    m_input_segments[i] = &m_in_segments[i];
            }

            void apply_angle_regularization() {

                m_angle_regularizer_ptr = std::make_shared<Angle_regularizer>(m_input_segments, m_parameters);
                
                m_angle_regularizer_ptr->make_silent(m_silent);
                m_angle_regularizer_ptr->regularize();
            }

            void apply_ordinate_regularization() {

                Tree *tree_pointer = m_angle_regularizer_ptr->get_tree_pointer();
                m_ordinate_regularizer_ptr = std::make_shared<Ordinate_regularizer>(m_input_segments, tree_pointer, m_parameters);

                m_ordinate_regularizer_ptr->make_silent(m_silent);
                m_ordinate_regularizer_ptr->regularize();
            }

            template<typename SegmentRange, typename SegmentMap>
            void update_input_segments(SegmentRange &input_segments, const SegmentMap &segment_map) {

                using Segment_iterator = typename SegmentRange::iterator;
                assert(m_input_segments.size() == input_segments.size());

                size_t i = 0;
                for (Segment_iterator sit = input_segments.begin(); sit != input_segments.end(); ++sit, ++i)
                    put(segment_map, *sit, m_input_segments[i]->get());

                assert(i == m_input_segments.size());
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H