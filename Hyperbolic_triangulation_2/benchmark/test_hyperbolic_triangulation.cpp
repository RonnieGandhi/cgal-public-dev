#include <fstream>

#include <CGAL/Hyperbolic_random_points_in_disc_2.h>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_triangulation_traits_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Hyperbolic_triangulation_traits_2<K> Gt;

typedef K::Point_2 Point_2;
typedef K::FT FT;

typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt> HDt;

int main(int argc, char *argv[])
{  
  FT r = 1;
  FT eps = 0;
  for(int k = 0; k < 2; k++) {
    if(k == 0) {
      eps = 1e-3;
    }
    if(k == 1) {
      eps = 1e-7;
    }
    
    int trials_nb = 10;
    int start_nb = 10000;
    
    std::cout << std::endl << eps << std::endl << std::endl;
    
    for(int nb = start_nb, k = 0; k < 4; nb = nb*10, k++) {
      
      std::vector< std::vector<Point_2> > pts(trials_nb);
      for(int i = 0; i < trials_nb; i++) {
        if(argc > 1 && argv[1][0] == 'e') {
          Random_points_in_disc_2<K>(pts[i], nb, i, eps);
        } else {
          Hyperbolic_random_points_in_disc_2<K>(pts[i], nb, i, eps);
        }
      }
      
      
      double average_time = 0;
      double average_nb = 0;
      double average_nb_of_edges = 0;
      for(int trials = 0; trials < trials_nb; trials++) {
        
        HDt hdt = HDt(Gt(r));
        
        CGAL::Timer timer;
        timer.start();
        
        hdt.insert(pts[trials].begin(), pts[trials].end());
        
        timer.stop();
        
        average_time += timer.time();
        timer.reset();
        
        average_nb += hdt.number_of_vertices();
        
        average_nb_of_edges = 0;
        for(HDt::Finite_edges_iterator eit = hdt.finite_edges_begin(); eit != hdt.finite_edges_end(); ++eit) {
          average_nb_of_edges++;
        }
      }
      average_time = average_time/trials_nb;
      average_nb = average_nb/trials_nb;
      
      std::cout << "H^2" << std::endl;
      std::cout << "Radius: " << r << std::endl;
      std::cout << "Eps: " << eps << std::endl;
      std::cout << "Number of points: " << average_nb << std::endl;
      std::cout << "Time: " << average_time << std::endl;
      std::cout << "number of edges " << average_nb_of_edges << std::endl;
      std::cout << std::endl;
    } 
  }
  
  return 0;
}
