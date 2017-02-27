

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


typedef CORE::Expr                                                                  NT;
typedef CGAL::Cartesian<NT>                                                         Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>         Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>                Triangulation;
typedef Triangulation::Face_iterator                                                Face_iterator;
typedef Triangulation::Vertex_handle 												Vertex_handle;
typedef Triangulation::Point 														Point;
typedef Traits::Side_of_fundamental_octagon                                         Side_of_fundamental_octagon;

typedef CGAL::Cartesian<double>::Point_2                                            Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                              Creator;

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number_of_iterations]" << endl;
        return -1;
    }

    int iter = atoi(argv[1]);
    Side_of_fundamental_octagon pred;  

    int min = 1000;
    int max = -1;
    double mean = 0.0;
    for (int j = 0; j < iter; j++) {
        cout << "Iteration " << (j+1) << "/" << iter << "..." << endl;
        
        vector<Point_double> v;
        Hyperbolic_random_points_in_disc_2_double(v, 500, -1, 0.159);

        Triangulation tr;  
        tr.insert_dummy_points(true);
        assert(tr.is_valid(true));
        
        int cnt = 0;
        int idx = 0;
        do {
            Point_double pt = v[idx++];
            if (pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
                tr.insert(Point(pt.x(), pt.y()));
                cnt++;
            }
        } while (tr.number_of_dummy_points() > 0);
        assert(tr.is_valid());
        
        if (cnt > max)
            max = cnt;
        if (cnt < min)
            min = cnt;
        mean += cnt;
    }
    mean /= (double)iter;

    cout << "Finished " << iter << " iterations!" << endl;
    cout << "Minimum number of points inserted: " << min << endl;
    cout << "Maximum number of points inserted: " << max << endl;
    cout << "Average number of points inserted: " << mean << endl;
    return 0;
}