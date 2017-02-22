

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_random_points_in_disc_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

#include <CGAL/Timer.h>



typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Hyperbolic_octagon_translation_matrix<NT>                               Octagon_matrix;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_fundamental_octagon                                     Side_of_fundamental_octagon;

typedef CGAL::Cartesian<double>::Point_2                                        Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                          Creator;


int main(void) {    

    Side_of_fundamental_octagon pred;

    int N = 10000;

    std::vector<Point> pts;
    cout << "---- for best results, make sure that you have compiled me in Release mode ----" << endl;
    cout << "generating " << N << " random points in the octagon..." << endl;    
    vector<Point_double> v;
    Hyperbolic_random_points_in_disc_2_double(v, 40*N, -1);

    int cnt = 0;
    int idx = 0;
    do {    
        Point pt = Point(v[idx].x(), v[idx].y());
        if (pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
            pts.push_back(pt);
            cnt++;
        } 
        idx++;
    } while (cnt < N && idx < v.size());

    if (cnt < N) {
        cout << "FAILED! Exiting..." << endl;
        return -1;
    }

    cout << "DONE!" << endl;

    cout << "inserting into triangulation with rational dummy points..." << endl;
    Triangulation tr2;
    tr2.insert_dummy_points(true);  

    CGAL::Timer t2;
    t2.start();
    for (int j = 0; j < pts.size(); j++) {
        //cout << "   now at " << j << endl;
        Vertex_handle vh = tr2.insert(pts[j]);
    }
    t2.stop();
    cout << " DONE!" << endl;
            
    cout << "inserting into triangulation with exact dummy points..." << endl;
    Triangulation tr1;
    tr1.insert_dummy_points(false);  

    CGAL::Timer t1;
    t1.start();
    for (int j = 0; j < pts.size(); j++) {
        cout << "   now at " << j << endl;
        Vertex_handle vh = tr1.insert(pts[j]);
    }
    t1.stop();
    cout << " DONE!" << endl;

    cout << "Triangulation with exact    dummy points: #vertices = " << tr1.number_of_vertices() << ", insertion time = " << t1.time() << endl;
    cout << "Triangulation with rational dummy points: #vertices = " << tr2.number_of_vertices() << ", insertion time = " << t2.time() << endl;

    return 0;
}
