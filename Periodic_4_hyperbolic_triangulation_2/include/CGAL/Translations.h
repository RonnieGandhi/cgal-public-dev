#ifndef CGAL_HYPERBOLIC_TRANSLATIONS_2_H
#define CGAL_HYPERBOLIC_TRANSLATIONS_2_H

#include <CGAL/Hyperbolic_isometry_2.h>

template<typename Translation>
struct Element
{
    typedef typename std::list<Element> List;
    typedef CGAL::Circulator_from_container<List> Circulator;
    
    Translation g;
    
    // circulator iterator to an inverse translation in the list
    Circulator inverse;
};

template<typename Gt>
class Translations
{
public:
    typedef typename Gt::FT FT;
    typedef typename std::complex<FT> complex;
    typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
    
    typedef Element<Hyperbolic_isometry> Element_t;
    typedef std::list<Element_t> List;
    typedef typename List::iterator List_iterator;
    typedef CGAL::Circulator_from_container<List> Circulator;
    
    typedef std::pair<Hyperbolic_isometry, int> Node;
    typedef std::vector<Node> Vector;
    typedef typename Vector::iterator Vector_iterator;
    
    
    static Hyperbolic_isometry& a()
    {
        compute();
        return g[0].first;
    }
    
    static Hyperbolic_isometry& b()
    {
        compute();
        //return g[1].first; // This was the original
        return g[5].first;
    } 
    
    static Hyperbolic_isometry& c()
    {
        compute();
        //return g[5].first; // This was the initial one
        return g[2].first;
    }
    
    static Hyperbolic_isometry& d()
    {
        compute();
        //return g[6].first; // This was the initial one
        return g[7].first;
    }
    
    static const Vector& get_vector_of_translations()
    {
        compute();
        return g;
    }
    
    static Vector_iterator vector_begin()
    {
        compute();
        return g.begin();
    }
    
    static Vector_iterator vector_end()
    {
        compute();
        return g.end();
    }
    
    static List_iterator list_begin()
    {
        compute();
        return l.begin();
    }
    
    static List_iterator list_end()
    {
        compute();
        return l.end();
    }
    
    static List& list()
    {
        compute();
        return l;
    }
    
private:
    
    static void compute_g()
    {
        // This is the x-coordinate (or real part) of the point on the unit circle where the line
        // with angle pi/8 (wrt the positive real half-axis) meets it.
        
        const FT k1 = (FT(2) + CGAL::sqrt(2.))/FT(2);
        const FT k2 = CGAL::sqrt(CGAL::sqrt(2.));
        const FT k3 = (CGAL::sqrt(2.)*k2)/FT(2);
        
        std::complex<FT> m(k1, k1);
        std::complex<FT> n(k2*k1, -k3);

        std::cout << "m = " <<  m << std::endl;
        std::cout << "n = " <<  n << std::endl;
        
        g.resize(8);
        
        // So, the following lines defining the group members don't seem to be correct -- it seems like they identify not
        // opposite sides, but instead sides which are on both sides of another one. Refer to page 91 in Manuel's thesis.
        // This is the case presented to the right of Figure 4.7, we want the left side.
        
        /*
         
         // a
         g[0].first = Hyperbolic_isometry(conj(m), conj(n));
         g[0].second = 2;
         
         // b
         g[3].first = Hyperbolic_isometry(m, -n);
         g[3].second = 1;
         
         // c
         g[4].first = Hyperbolic_isometry(conj(m), -conj(n));
         g[4].second = 6;
         
         // d
         g[7].first = Hyperbolic_isometry(m, n);
         g[7].second = 5;
         
         
         
         int index = g[0].second;
         g[index].first = g[0].first.inverse();
         g[index].second = 0;
         
         index = g[3].second;
         g[index].first = g[3].first.inverse();
         g[index].second = 3;
         
         index = g[4].second;
         g[index].first = g[4].first.inverse();
         g[index].second = 4;
         
         index = g[7].second;
         g[index].first = g[7].first.inverse();
         g[index].second = 7;
         */
        
        // ----------------------------------------------------------//
        
        // Here begins my own implementation -- let's see what happens.
        
        // a
         g[0].first = Hyperbolic_isometry(conj(m), conj(n));
         g[0].second = 2;
         
         // b
         g[3].first = Hyperbolic_isometry(m, -n);
         g[3].second = 1;
         
         // c
         g[4].first = Hyperbolic_isometry(conj(m), -conj(n));
         g[4].second = 6;
         
         // d
         g[7].first = Hyperbolic_isometry(m, n);
         g[7].second = 5;
        
    }
    
    
    static void compute_l()
    {
        l.resize(g.size());
        
        std::vector<Circulator> aux_list;
        aux_list.reserve(8);
        
        for(List_iterator li = l.begin(); li != l.end(); li++) {
            aux_list.push_back( Circulator(&l, li) );
        }
        
        for(typename List::size_type i = 0; i < aux_list.size(); i++) {
            aux_list[i]->g = g[i].first;
            aux_list[i]->inverse = aux_list[g[i].second];
        }
    }
    
    static void compute()
    {
        static bool computed = false;
        if(!computed) {
            compute_g();
            compute_l();
            computed = true;
        }
    }
    
    static Vector g;
    
    static List l;
};

// default initialization

template<typename Gt>
typename Translations<Gt>::Vector
Translations<Gt>::g;

template<typename Gt>
typename Translations<Gt>::List
Translations<Gt>::l;

#endif // CGAL_HYPERBOLIC_TRANSLATIONS_2_H
