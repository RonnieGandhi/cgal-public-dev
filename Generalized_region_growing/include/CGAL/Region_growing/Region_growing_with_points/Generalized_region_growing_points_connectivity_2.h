#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H

#include <CGAL/Iterator_range.h>

template <class Traits>
class Generalized_region_growing_points_connectivity_2 {
public:
    typedef typename Traits::Element Element;
    typedef typename Traits::Input_range Input_range;
    typedef typename Traits::Element_map Element_map;
    typedef typename Traits::Normal_map Normal_map;

    typedef typename Input_range::iterator Iterator;
    typedef typename std::vector<Iterator> Neighbors;
    typedef typename Iterator_range<Neighbors::iterator> Neighbor_range;
    
    Generalized_region_growing_points_connectivity_2(Input_range input_range) :
        m_input_range(input_range) {}

    void get_neighbors(const Iterator& elem_iter, Neighbor_range& output) const {
        m_neighbors.clear();
        Element elem1 = get(elem_map, elem_iter);
        /*
        for (Iterator iter = m_input_range.begin(); iter != m_input_range.end(); ++iter) {
            Element elem2 = get(elem_map, iter);
            ...
            m_neighbors.push_back(iter);
        }
        */
       output = Neighbor_range(m_neighbors.begin(), m_neighbors.end());
    }
private:
    Neighbors m_neighbors;
    Element_map m_elem_map;
    Normal_map m_normal_map;
    Input_range m_input_range;
}

#endif
