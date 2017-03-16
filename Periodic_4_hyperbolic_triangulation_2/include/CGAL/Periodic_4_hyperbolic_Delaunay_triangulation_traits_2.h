// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:  $
//
//
// Author(s)     : 	Iordan Iordanov
// 					Monique Teillaud
//


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
#include <CGAL/exact_complex.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"


namespace CGAL {


template < class K, class Predicate_ >
class Hyperbolic_traits_with_offsets_2_adaptor
{
	typedef K 				Kernel;
	typedef Predicate_ 		Predicate;

	typedef typename Kernel::Point_2    	Point;
	typedef typename Kernel::FT 			FT;
	typedef typename Kernel::Offset  	 	Offset;

	// Use the construct_point_2 predicate from the kernel to convert the periodic points to Euclidean points
	typedef typename Kernel::Construct_point_2        			Construct_point_2;

public:
	typedef typename Predicate::result_type           			result_type;


	Hyperbolic_traits_with_offsets_2_adaptor() { }

	result_type operator()(	const Point& p0, 	const Point& p1,
							const Offset& o0, 	const Offset& o1) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1));
	}
	result_type operator()(	const Point& p0, 	const Point& p1, 	const Point& p2,
							const Offset& o0, 	const Offset& o1, 	const Offset& o2) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2, 	const Point& p3,
							const Offset& o0, 	const Offset& o1,
							const Offset& o2, 	const Offset& o3) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3));
	}

	result_type operator()(	const Point& p0, 	const Point& p1) const
	{
		return Predicate()(p0, p1);
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2) 	const
	{
		return Predicate()(p0, p1, p2);
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2, 	const Point& p3) const
	{
		return Predicate()(p0, p1, p2, p3);
	}

private:
	Point pp(const Point &p, const Offset &o) const
	{
		return o.apply(p);
	}

};



template < typename K, typename Construct_point_base>
class Periodic_4_hyperbolic_construct_point_2 : public Construct_point_base
{
	typedef K Kernel;

public:
	typedef typename Kernel::Point_2         Point;
	typedef typename Kernel::Offset          Offset;

	typedef Point                            result_type;

	Periodic_4_hyperbolic_construct_point_2() { }

	Point operator() ( const Point& p, const Offset& o ) const
	{
		return o.apply(p);
	}

};



template <class R>
class Simple_circular_arc_2 {
	typedef typename R::FT 			FT;
	typedef typename R::Point_2 	Point;
	typedef typename R::Circle_2 	Circle;

private:
	Circle _c;
	Point _s, _t;

public:
	Simple_circular_arc_2() :
		_c(Point(FT(0),FT(0)), FT(0)), _s(FT(0),FT(0)), _t(FT(0),FT(0)) {}
	
	Simple_circular_arc_2(Circle c, Point source, Point target) :
		_c(c), _s(source), _t(target) {}

	Circle circle() const {
		return _c;
	}

	Point source() const {
		return _s;
	}

	Point target() const {
		return _t;
	}

	FT squared_radius() const {
		return _c.squared_radius();
	}

	Point center() const {
		return _c.center();
	}

	Bbox_2 bbox(void) const {
    	return typename R::Construct_bbox_2()(*this);
  	}

};



template< class R >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 : public R {

typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<R>  Self;  

public:

	typedef typename R::FT          								FT;
	typedef Hyperbolic_octagon_word_4<unsigned short int, FT>		Offset;
	typedef typename R::Point_2     								Point_2;
	typedef Point_2                 								Point;
	typedef typename R::Circle_2    								Circle_2;
	typedef typename R::Line_2      								Euclidean_line_2;
	typedef boost::variant<Circle_2,Euclidean_line_2>    			Euclidean_circle_or_line_2; 
	typedef Simple_circular_arc_2<R>         						Circular_arc_2;
	typedef typename R::Segment_2                       			Euclidean_segment_2; //only used internally here
	typedef boost::variant<Circular_arc_2, Euclidean_segment_2>  	Hyperbolic_segment_2;

	// Wrappers for the offset adapter
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Compare_x_2>                 Compare_x_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Compare_y_2>                 Compare_y_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Orientation_2>               Orientation_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Compare_distance_2>          Compare_distance_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Side_of_oriented_circle_2>   Side_of_oriented_circle_2;

	
	

	// only kept for demo to please T2graphicsitems
	typedef Euclidean_segment_2  							Line_segment_2;
	typedef Hyperbolic_segment_2 							Segment_2;

	// the following types are only used internally in this traits class, 
	// so they need not be documented, and they don't need _object()
	typedef typename R::Collinear_2                			Euclidean_collinear_2;
	typedef typename R::Construct_bisector_2       			Construct_Euclidean_bisector_2;
	typedef typename R::Construct_midpoint_2       			Construct_Euclidean_midpoint_2;
	typedef typename R::Compute_squared_distance_2 			Compute_squared_Euclidean_distance_2;
	typedef typename R::Has_on_bounded_side_2 				Has_on_bounded_side_2;

	typedef typename R::Less_x_2                   			Less_x_2;
	typedef typename R::Less_y_2                   			Less_y_2;
			
public:

	class Construct_hyperbolic_segment_2 {
		typedef exact_complex<FT> 	cplx;
	public:
		typedef Segment_2 result_type;

		Construct_hyperbolic_segment_2() {}

		Segment_2 operator()(const Point_2& p1, const Point_2& p2) const {
			cplx p(p1), q(p2);
			cplx O(0,0);
			cplx inv;
			if (p == O) {
				inv = q.invert_in_unit_circle();
			} else {
				inv = p.invert_in_unit_circle();
			}

			Point ip(inv.real(), inv.imag());

			if (Orientation_2()(p1, p2, ip) == COLLINEAR) {
				Euclidean_segment_2 seg(p1, p2);
				return seg;
			} else {
				Circle_2 c(p1, p2, ip);
				if(Orientation_2()(p1, p2, c.center()) == LEFT_TURN) {
					return Circular_arc_2(c, p1, p2);
				}
				return Circular_arc_2(c, p2, p1);
			}

		}
		
	}; // end Construct_hyperbolic_segment_2
	
	Construct_hyperbolic_segment_2
		construct_hyperbolic_segment_2_object() const
	{ return Construct_hyperbolic_segment_2(); }


	// wrong names kept for demo
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, Construct_hyperbolic_segment_2> Construct_segment_2;
	Construct_segment_2
	construct_segment_2_object() const
	{ return Construct_segment_2(); }
	
	
	Compare_x_2 
		compare_x_2_object() const 
	{ return Compare_x_2();} 
	
	Compare_y_2 
		compare_y_2_object() const 
	{ return Compare_y_2();} 
	
	Orientation_2
		orientation_2_object() const
	{ return Orientation_2();}
	
	Side_of_oriented_circle_2
		side_of_oriented_circle_2_object() const
	{ return Side_of_oriented_circle_2(); }
	


	class Construct_hyperbolic_circle_2 {

		typedef exact_complex<FT> 	cplx;

	public: 
		Construct_hyperbolic_circle_2() {}

		Circle_2 operator()(Point_2 hcenter, Point_2 p) {
			
			Point o(0,0);

			if (hcenter == o) {
				return Circle_2(o, squared_distance(o, p));
			} else if (Orientation_2()(hcenter, p, o) != COLLINEAR) {
				Euclidean_line_2 ell(hcenter, o);
				
				cplx p1(hcenter), p2(p);
				cplx inv;
				if (p1 == cplx(0,0)) {
					inv = p2.invert_in_unit_circle();
				} else {
					inv = p1.invert_in_unit_circle();
				}
				Point ip(inv.real(), inv.imag());

				Circle_2 schl(hcenter, p, ip);
				Euclidean_line_2 line_through_p(schl.center(), p);
				Euclidean_line_2 tangent_at_p = line_through_p.perpendicular(p);

				// assume that ell := ax + by + c = 0 and tangent_at_p := dx + ey + f = 0
				FT a = ell.a(), b = ell.b(), c = ell.c();
				FT d = tangent_at_p.a(), e = tangent_at_p.b(), f = tangent_at_p.c();
                FT px, py;
                if (a == FT(0)) {
                	px = (e*c - b*f)/(d*b);
                	py = -c/b;
                } else if (b == FT(0)) {
                	px = -c/a;
                	py = (c*d - a*f)/(a*e);
                } else if (d == FT(0)) {
                	px = (b*f - e*c)/(e*a);
                	py = -f/e;
                } else if (e == FT(0)) {
                	px = -f/d;
                	py = (a*f - d*c)/(d*b);
                } else {
                	py = (c*d - a*f)/(a*e - d*b);	
                	px = (-c -b*py)/a;
                }
				Point intersection(px, py);
				return Circle_2(intersection, squared_distance(intersection, p));
			} else {  // if the given points and the origin are collinear, we need to treat them differently
				cplx hcinv = cplx(hcenter).invert_in_unit_circle();
				Point ip(hcinv.real(), hcinv.imag());
				Point mp = midpoint(hcenter, ip);
				Circle_2 tmpc(mp, hcenter);
				cplx res = cplx(p).invert_in_circle(tmpc);
				Point pres(res.real(), res.imag());
				return Circle_2(p, pres);
			}
		}
	};




	class Construct_inexact_hyperbolic_bisector_2 {

		typedef exact_complex<FT> 	cplx;

	public:      
		Construct_inexact_hyperbolic_bisector_2() 
			{}
		
		Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q) const
		{

			Origin o; 
			Point_2 po = Point_2(o);
			double np = sqrt(to_double(p.x()*p.x() + p.y()*p.y()));
			double nq = sqrt(to_double(q.x()*q.x() + q.y()*q.y()));
			if ( fabs(np - nq) < 1e-16 ) {      
				Euclidean_line_2 ebl = Construct_Euclidean_bisector_2()(p, q);
				pair<Point_2, Point_2> pts = Construct_inexact_intersection_2()(ebl, Circle_2(Point_2(0, 0), 1));
				return Euclidean_segment_2(pts.first, pts.second);
			}

			Circle_2 c1 = Construct_hyperbolic_circle_2()(p, q);
			Circle_2 c2 = Construct_hyperbolic_circle_2()(q, p);

			pair<Point_2, Point_2> res = Construct_inexact_intersection_2()(c1, c2);
			Point_2 p1 = res.first, p2 = res.second;

			cplx cp1(p1), cp2(p2);
			cplx inv;
			if (cp1 == cplx(0,0)) {
				inv = cp2.invert_in_unit_circle();
			} else {
				inv = cp1.invert_in_unit_circle();
			}
			
			Point cpi(inv.real(), inv.imag());

			Circle_2 c(p1, p2, cpi);

			res = Construct_inexact_intersection_2()(c, Circle_2(Point(0,0), 1));
			p1 = res.first;
			p2 = res.second;

			if(Orientation_2()(p1, p2, c.center()) == LEFT_TURN) {
				return Circular_arc_2(c, p1, p2);
			}
			return Circular_arc_2(c, p2, p1);

		}
	}; // end Construct_hyperbolic_bisector_2
	
	Construct_inexact_hyperbolic_bisector_2
	construct_inexact_hyperbolic_bisector_2_object() const
	{ return Construct_inexact_hyperbolic_bisector_2(); }





	class Construct_hyperbolic_bisector_2 {

		typedef exact_complex<FT> 	cplx;

	public:      
		Construct_hyperbolic_bisector_2() 
			{}
		
		Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q) const
		{

			Origin o; 
			Point_2 po = Point_2(o);
			if ( Compare_distance_2()(po, p, q) == EQUAL ){      
				Euclidean_line_2 ebl = Construct_Euclidean_bisector_2()(p, q);
				pair<Point_2, Point_2> pts = Construct_intersection_2()(ebl, Circle_2(Point_2(0, 0), 1));
				return Euclidean_segment_2(pts.first, pts.second);
			}

			Circle_2 c1 = Construct_hyperbolic_circle_2()(p, q);
			Circle_2 c2 = Construct_hyperbolic_circle_2()(q, p);

			pair<Point_2, Point_2> res = Construct_intersection_2()(c1, c2);
			Point_2 p1 = res.first, p2 = res.second;

			cplx cp1(p1), cp2(p2);
			cplx inv;
			if (cp1 == cplx(0,0)) {
				inv = cp2.invert_in_unit_circle();
			} else {
				inv = cp1.invert_in_unit_circle();
			}
			
			Point cpi(inv.real(), inv.imag());

			Circle_2 c(p1, p2, cpi);

			res = Construct_intersection_2()(c, Circle_2(Point(0,0), 1));
			p1 = res.first;
			p2 = res.second;

			if(Orientation_2()(p1, p2, c.center()) == LEFT_TURN) {
				return Circular_arc_2(c, p1, p2);
			}
			return Circular_arc_2(c, p2, p1);

		}
	}; // end Construct_hyperbolic_bisector_2
	
	Construct_hyperbolic_bisector_2
	construct_hyperbolic_bisector_2_object() const
	{ return Construct_hyperbolic_bisector_2(); }
	
	Construct_Euclidean_bisector_2
	construct_Euclidean_bisector_2_object() const
	{ return Construct_Euclidean_bisector_2(); }	



	class Construct_intersection_2 {
	public:
		Construct_intersection_2() {}

		Point_2 operator()(Euclidean_line_2 ell1, Euclidean_line_2 ell2) {
			//cout << "intersecting lines..." << endl;
			if (ell1.b() == FT(0)) {
				swap(ell1, ell2);
			}
			
			CGAL_assertion(ell1.b() != FT(0));
			if (ell2.b() != FT(0)) {
				CGAL_assertion( ell1.a()/ell1.b() != ell2.a()/ell2.b() );
			}

			FT lambda1 = -ell1.a()/ell1.b();
			FT mu1     = -ell1.c()/ell1.b();
			FT x = ( -ell2.c() - mu1*ell2.b() )/( ell2.a() + lambda1*ell2.b() );
			FT y = lambda1*x + mu1;
			return Point_2(x, y);
		}

		std::pair<Point_2, Point_2> operator()(Euclidean_line_2 ell, Circle_2 c) {
			//cout << "intersecting line and circle: " << ell.a() << "x + " << ell.b() << "y + " << ell.c() << " = 0, c = " << c.center() << ", r_sq = " << c.squared_radius() << endl;
			if (ell.b() == FT(0)) {
			//	cout << "  degenerate case!" << endl;
				FT p 	= c.center().x();
				FT q 	= c.center().y(); 
				FT y1  	= q + CGAL::sqrt(c.squared_radius() - p*p);
				FT y2  	= q - CGAL::sqrt(c.squared_radius() - p*p);
				Point_2 p1(FT(0), y1);
				Point_2 p2(FT(0), y2);
				return make_pair(p1, p2);
			}

			//cout << "  non-degenerate case" << endl;
			FT lambda = -ell.a()/ell.b();
			FT mu 	  = -ell.c()/ell.b();
			FT p 	  = c.center().x();
			FT q 	  = c.center().y(); 
			FT A = FT(1) + lambda*lambda;
			FT B = FT(2)*( lambda * mu - lambda*q - p);
			FT C = p*p + mu*mu + q*q - c.squared_radius() - FT(2)*q*mu;
			FT Delta = B*B - FT(4)*A*C;
			//cout << "  computed quantities" << endl;
			//cout << "  A     = " << A << endl;
			//cout << "  B     = " << B << endl;
			//cout << "  C     = " << C << endl;
			//cout << "  Delta = " << Delta << endl;
			//CGAL_assertion(Delta >= FT(0));
			// if (Delta == FT(0)) {
			// 	cout << "  one double root" << endl;
			// 	FT x = -B/(FT(2)*A);
			// 	FT y = lambda*x + mu;
			// 	Point_2 sol(x, y);
			// 	return make_pair(sol, sol);
			// } else {
				FT x1 = (-B + CGAL::sqrt(Delta))/(FT(2)*A);
				//cout << "    x1 = " << x1 << endl;
				FT x2 = (-B - CGAL::sqrt(Delta))/(FT(2)*A);
				//cout << "    x2 = " << x2 << endl;
				FT y1 = lambda*x1 + mu;
				//cout << "    y1 = " << y1 << endl;
				FT y2 = lambda*x2 + mu;
				//cout << "    y2 = " << y2 << endl;
				Point_2 sol1(x1, y1);
				Point_2 sol2(x2, y2);
				return make_pair(sol1, sol2);
			//}
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c, Euclidean_line_2 ell) {
			return operator()(ell, c);
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c1, Circle_2 c2) {
			//cout << "intersecting circles: c1 = " << c1.center() << ", r1_sq = " << c1.squared_radius() << ", c2 = " << c2.center() << ", r2_sq = " << c2.squared_radius() << endl;
			
			FT xa = c1.center().x(), ya = c1.center().y();
			FT xb = c2.center().x(), yb = c2.center().y();
			//cout << "   centers" << endl;
			FT d2 = squared_distance(c1.center(), c2.center());
			//cout << "   centers squared distance" << endl;
			FT ra = CGAL::sqrt(c1.squared_radius());
			//cout << "   radius 1" << endl;
			FT rb = CGAL::sqrt(c2.squared_radius());
			//cout << "   radius 2" << endl;
			FT K  = CGAL::sqrt(((ra+rb)*(ra+rb)-d2)*(d2-(ra-rb)*(ra-rb)))/FT(4); 
			//cout << "   quantities ready" << endl;

			FT xbase = (xb + xa)/FT(2) + (xb - xa)*(ra*ra - rb*rb)/d2/FT(2);
			FT xdiff = FT(2)*(yb - ya)*K/d2;
			FT x1 = xbase + xdiff;
			FT x2 = xbase - xdiff;
			//cout << "   x computed" << endl;

			FT ybase = (yb + ya)/FT(2) + (yb - ya)*(ra*ra - rb*rb)/d2/FT(2);
			FT ydiff = FT(-2)*(xb - xa)*K/d2;
			FT y1 = ybase + ydiff;
			FT y2 = ybase - ydiff;
			//cout << "   y computed" << endl;

			Point_2 res1(x1, y1);
			Point_2 res2(x2, y2);
			return make_pair(res1, res2);

			// FT p1 = c1.center().x();
			// FT p2 = c2.center().x();
			// FT q1 = c1.center().y();
			// FT q2 = c2.center().y();

			// FT a = FT(2)*(p2 - p1);
			// FT b = FT(2)*(q2 - q1);
			// FT c = p1*p1 + q1*q1 + c2.squared_radius() - p2*p2 - q2*q2 - c1.squared_radius();
			// Euclidean_line_2 ell(a, b, c);
			// return operator()(ell, c1);
		}


		Point_2 operator()(Hyperbolic_segment_2 s1, Hyperbolic_segment_2 s2) {
			if (Circular_arc_2* c1 = boost::get<Circular_arc_2>(&s1)) {
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(c1->circle(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					pair<Point_2, Point_2> res = operator()(c1->circle(), ell2->supporting_line());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				}
			} else {
				Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(ell1->supporting_line(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;	
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					Point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
					CGAL_assertion(CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1));
					return p1;
				}
			}
		}

	};


	Construct_intersection_2
	construct_intersection_2_object() const {
		return Construct_intersection_2();
	}




	class Construct_inexact_intersection_2 {
	public:
		Construct_inexact_intersection_2() {}

		Point_2 operator()(Euclidean_line_2 ell1, Euclidean_line_2 ell2) {
			//cout << "inexactly intersecting lines..." << endl;

			if (fabs(to_double(ell1.b())) < 1e-16) {
				swap(ell1, ell2);
			}
			
			double a1 = to_double(ell1.a()), b1 = to_double(ell1.b()), c1 = to_double(ell1.c());
			double a2 = to_double(ell2.a()), b2 = to_double(ell2.b()), c2 = to_double(ell2.c());

			CGAL_assertion(fabs(b1) > 1e-16);
			if (fabs(b2) > 1e-16) {
				CGAL_assertion( fabs(a1/b1 - a2/b2) > 1e-16 );
			}

			double lambda1 = -a1/b1;
			double mu1     = -c1/b1;
			double x = ( -c2 - mu1*b2 )/( a2 + lambda1*b2 );
			double y = lambda1*x + mu1;
			return Point_2(x, y);
		}

		std::pair<Point_2, Point_2> operator()(Euclidean_line_2 ell, Circle_2 cc) {
			double a = to_double(ell.a()), b = to_double(ell.b()), c = to_double(ell.c());
			double p = to_double(cc.center().x()), q = to_double(cc.center().y()), r2 = to_double(cc.squared_radius());
			//cout << "inexactly intersecting line and circle: " << a << "x + " << b << "y + " << c << " = 0, c = (" << p << "," << q << "), r_sq = " << r2 << endl;
			
			double A, B, C, D;
			double x1, y1, x2, y2;
			if (fabs(a) < 1e-16) {
				y1 = -c/b;  y2 = -c/b;
				A = b*p;
				D = -b*b*q*q + b*b*r2 - 2.*b*c*q - c*c;
				//cout << "A     = " << A 	<< endl;
				//cout << "D     = " << D  	<< endl;
				x1 = (A + sqrt(D))/b;
				x2 = (A - sqrt(D))/b;
			} else if (fabs(b) < 1e-16) {
				x1 = -c/a;  x2 = -c/a;
				A = q*a;
				D = -a*a*p*p + r2*a*a - 2.*a*c*p - c*c;
				//cout << "A     = " << A 	<< endl;
				//cout << "D     = " << D  	<< endl;
				y1 = (A + sqrt(D))/a;
				y2 = (A - sqrt(D))/a;
			} else {
				A = a*a*q - a*b*p-b*c;
				C = (-b*q - c)*a*a + b*b*p*a;
				D = -a*a*( b*b*q*q + 2.*q*(p*a + c)*b - b*b*r2 + (p*p - r2)*a*a + 2.*a*c*p + c*c );
				B = a*a + b*b;

				//cout << "A     = " << A 	<< endl;
				//cout << "B     = " << B 	<< endl;
				//cout << "C     = " << C 	<< endl;
				//cout << "D     = " << D  	<< endl;

				y1 = (A + sqrt(D))/B;
				y2 = (A - sqrt(D))/B;
				x1 = (C - b*sqrt(D))/(a*(a*a + b*b));
				x2 = (C + b*sqrt(D))/(a*(a*a + b*b));
			}

			Point_2 p1(x1, y1);
			Point_2 p2(x2, y2);

			return make_pair(p1, p2);
			// double A, B, C, Delta, lambda, mu;
			// if (fabs(b) < 1e-16) {
			// 	cout << "  degenerate case!" << endl;
			// 	double mu = -c/a;
			// 	A = 1.;
			// 	B = -2.*q;
			// 	C = q*q - r2 + (mu-p)*(mu-p);
				
			// 	Delta = B*B - 4.*A*C;
			// 	cout << "mu    = " << mu 	<< endl;
			// 	cout << "A     = " << A 	<< endl;
			// 	cout << "B     = " << B 	<< endl;
			// 	cout << "C     = " << C 	<< endl;
			// 	cout << "Delta = " << Delta << endl;

			// 	double y1 = (-B + sqrt(Delta))/(2.*A);
			// 	double y2 = (-B - sqrt(Delta))/(2.*A);
			// 	double x = mu;

			// 	cout << "    x1 = " << x << endl;
			// 	cout << "    y1 = " << y1 << endl;
			// 	cout << "    x2 = " << x << endl;
			// 	cout << "    y2 = " << y2 << endl;

			// 	Point_2 sol1(x, y1);
			// 	Point_2 sol2(x, y2);
			// 	return make_pair(sol1, sol2);
			// } else {
			// 	cout << "  non-degenerate case" << endl;
			// 	double lambda = -a/b;
			// 	double mu 	  = -c/b;
			// 	A = 1. + lambda*lambda;
			// 	B = 2.*( lambda*mu - lambda*p - q);
			// 	C = mu*mu + p*p + q*q - r2 - 2.*mu*p;
			
			// 	Delta = B*B - 4.*A*C;
					
			// 	cout << "A     = " << A 	<< endl;
			// 	cout << "B     = " << B 	<< endl;
			// 	cout << "C     = " << C 	<< endl;
			// 	cout << "Delta = " << Delta << endl;

			// 	double y1 = (-B + sqrt(Delta))/(2.*A);
			// 	double y2 = (-B - sqrt(Delta))/(2.*A);
			// 	double x1 = lambda*y1 + mu;
			// 	double x2 = lambda*y2 + mu;

			// 	cout << "    x1 = " << x1 << endl;
			// 	cout << "    y1 = " << y1 << endl;
			// 	cout << "    x2 = " << x2 << endl;
			// 	cout << "    y2 = " << y2 << endl;

			// 	Point_2 sol1(x1, y1);
			// 	Point_2 sol2(x2, y2);
			// 	return make_pair(sol1, sol2);
			
			// }
			
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c, Euclidean_line_2 ell) {
			return operator()(ell, c);
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c1, Circle_2 c2) {
			//cout << "inexactly intersecting circles..." << endl;
			double xa = to_double(c1.center().x()), ya = to_double(c1.center().y());
			double xb = to_double(c2.center().x()), yb = to_double(c2.center().y());
			//cout << "   centers" << endl;
			double d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
			//cout << "   centers squared distance" << endl;
			double ra = sqrt(to_double(c1.squared_radius()));
			//cout << "   radius 1" << endl;
			double rb = sqrt(to_double(c2.squared_radius()));
			//cout << "   radius 2" << endl;
			double K  = sqrt(((ra+rb)*(ra+rb)-d2)*(d2-(ra-rb)*(ra-rb)))/4.; 
			//cout << "   quantities ready" << endl;

			double xbase = (xb + xa)/2. + (xb - xa)*(ra*ra - rb*rb)/d2/2.;
			double xdiff = 2.*(yb - ya)*K/d2;
			double x1 = xbase + xdiff;
			double x2 = xbase - xdiff;
			//cout << "   x computed" << endl;

			double ybase = (yb + ya)/2. + (yb - ya)*(ra*ra - rb*rb)/d2/2.;
			double ydiff = -2.*(xb - xa)*K/d2;
			double y1 = ybase + ydiff;
			double y2 = ybase - ydiff;
			//cout << "   y computed" << endl;

			Point_2 res1(x1, y1);
			Point_2 res2(x2, y2);
			//cout << "   returning!" << endl;
			return make_pair(res1, res2);
		}


		Point_2 operator()(Hyperbolic_segment_2 s1, Hyperbolic_segment_2 s2) {
			if (Circular_arc_2* c1 = boost::get<Circular_arc_2>(&s1)) {
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(c1->circle(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					pair<Point_2, Point_2> res = operator()(c1->circle(), ell2->supporting_line());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				}
			} else {
				Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(ell1->supporting_line(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;	
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					Point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
					CGAL_assertion(CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1));
					return p1;
				}
			}
		}

	};



	Construct_inexact_intersection_2
	construct_inexact_intersection_2_object() const {
		return Construct_inexact_intersection_2();
	}






	class Construct_hyperbolic_circumcenter_2_base {

	public:

		typedef Point_2 result_type;

		Construct_hyperbolic_circumcenter_2_base() {}

		Point_2 operator()(Point_2 p, Point_2 q, Point_2 r) {

			Hyperbolic_segment_2 s1 = Construct_hyperbolic_bisector_2()(p, q);
			Hyperbolic_segment_2 s2 = Construct_hyperbolic_bisector_2()(p, r);

			Point_2 rp = Construct_intersection_2()(s1, s2);
			return rp;
		}

	};


	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, Construct_hyperbolic_circumcenter_2_base> Construct_hyperbolic_circumcenter_2;


	Construct_hyperbolic_circumcenter_2
	construct_hyperbolic_circumcenter_2_object() const {
		return Construct_hyperbolic_circumcenter_2();
	}





	class Construct_inexact_hyperbolic_circumcenter_2_base {

	public:

		typedef Point_2 result_type;

		Construct_inexact_hyperbolic_circumcenter_2_base() {}

		Point_2 operator()(Point_2 p, Point_2 q, Point_2 r) {

			Hyperbolic_segment_2 s1 = Construct_inexact_hyperbolic_bisector_2()(p, q);
			Hyperbolic_segment_2 s2 = Construct_inexact_hyperbolic_bisector_2()(p, r);

			Point_2 rp = Construct_inexact_intersection_2()(s1, s2);
			return rp;
		}

	};


	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, Construct_inexact_hyperbolic_circumcenter_2_base> Construct_inexact_hyperbolic_circumcenter_2;


	Construct_inexact_hyperbolic_circumcenter_2
	construct_inexact_hyperbolic_circumcenter_2_object() const {
		return Construct_inexact_hyperbolic_circumcenter_2();
	}


	/****************************************************/
	class Side_of_hyperbolic_face_2 {
		
	public:
		typedef Bounded_side result_type;

		Side_of_hyperbolic_face_2() 
			{}



		template<class Face_handle, class Offset>
		Bounded_side operator()(const Point_2 p, Bounded_side sides[3], const Face_handle fh, const Offset o) const {

			Point_2 p1 = o.append(fh->offset(0)).apply(fh->vertex(0)->point());
			Point_2 p2 = o.append(fh->offset(1)).apply(fh->vertex(1)->point());
			Point_2 p3 = o.append(fh->offset(2)).apply(fh->vertex(2)->point());

			Bounded_side cp1 = side_of_segment_2(p,  p2, p3);
			sides[0] = cp1;
			if (cp1 == ON_BOUNDARY) {
				return ON_BOUNDARY;
			}

			Bounded_side cp2 = side_of_segment_2(p,  p3, p1);
			sides[1] = cp2;
			if (cp2 == ON_BOUNDARY) {
				return ON_BOUNDARY;
			}

			Bounded_side cp3 = side_of_segment_2(p,  p1, p2);
			sides[2] = cp3;
			if (cp3 == ON_BOUNDARY) {
				return ON_BOUNDARY;
			}			

			Bounded_side cs1 = side_of_segment_2(p1, p2, p3);
			Bounded_side cs2 = side_of_segment_2(p2, p3, p1);
			Bounded_side cs3 = side_of_segment_2(p3, p1, p2);


			// Cannot be on the boundary here.
			if (cs1 != cp1 || cs2 != cp2 || cs3 != cp3) {
				return ON_UNBOUNDED_SIDE;
			} else {
				return ON_BOUNDED_SIDE;	
			}


		}



		template<class Face_handle>
		Bounded_side operator()(const Point_2 p, Bounded_side sides[3], const Face_handle fh) const {
			return operator()(p, sides, fh, Offset());
		}


	private:
		Bounded_side side_of_segment_2(const Point_2 query, const Point_2 p, const Point_2 q) const {
			
			Point_2 o(0, 0);

			// Invert p or q through the unit circle.
			// The inversion depends on the distance from the origin, so to increase 
			// numerical stability we choose to invert the point further from (0,0).
			Point_2 inv;
			FT dp = squared_distance(o, p), dq = squared_distance(o, q);
			if (dq < dp) {
				inv = Point_2( p.x()/dp, p.y()/dp );
			} else {
				inv = Point_2( q.x()/dq, q.y()/dq );
			}

			// If iq is on the line defined by p and q, we need to work with a line instead of a circle
			if (orientation(p, q, inv) == COLLINEAR) {
				Orientation oquery = orientation(p, q, query);
				if (oquery == COLLINEAR) {
					return ON_BOUNDARY;
				} else if (oquery == LEFT_TURN) {
					return ON_BOUNDED_SIDE; 	// this is just a convention
				} else {
					return ON_UNBOUNDED_SIDE;
				}
			} else { // this means that we work in the circle
				Circle_2 c(p, q, inv);
				return c.bounded_side(query);
			}

		}

	};


	Side_of_hyperbolic_face_2
	side_of_hyperbolic_face_2_object() const {
		return Side_of_hyperbolic_face_2();
	}

	/****************************************************/




public:
		Periodic_4_hyperbolic_Delaunay_triangulation_traits_2() {}
	
		Periodic_4_hyperbolic_Delaunay_triangulation_traits_2(const Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 & other) {}
	
		Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 &operator=(const Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 &)
		{
			return *this;
		}



		class Side_of_fundamental_octagon {
		public:
			Side_of_fundamental_octagon() {}

			template <class Point_2_template>
			CGAL::Bounded_side operator()(Point_2_template p) {

				FT F2(2);
				FT qty = CGAL::sqrt(F2 + F2*CGAL::sqrt(F2));
				// The center of the Euclidean circle corresponding to the side s_1 (east)
				Point_2 CenterA ( qty/F2, FT(0) );
				Point_2 CenterB ( qty*CGAL::sqrt(F2)/FT(4), qty*CGAL::sqrt(F2)/FT(4) );

				// The squared radius of the Eucliden circle corresponding to the side s_1
				FT      Radius2 ( (CGAL::sqrt(F2) - FT(1)) / F2 );

				// Poincare disk (i.e., unit Euclidean disk)
				Circle_2 Poincare    ( Point(0, 0), 1*1 );

				// Euclidean circle corresponding to s_1
				Circle_2 EuclidCircA ( CenterA,     Radius2 );

				// Euclidean circle corresponding to s_2 (just rotate the center, radius is the same)
				Circle_2 EuclidCircBb( CenterB,     Radius2 );

				// This transformation brings the point in the first quadrant (positive x, positive y)
				FT x(FT(p.x()) > FT(0) ? p.x() : -p.x());
				FT y(FT(p.y()) > FT(0) ? p.y() : -p.y());

				// This brings the point in the first octant (positive x and y < x)
				if (y > x) {
					FT tmp = x;
					x = y;
					y = tmp;
				}

				// This tells us whether the point is in the side of the open boundary
				bool on_open_side = ( ( p.y() + tan(CGAL_PI / 8.) * p.x() ) < 0.0 );

				Point t(x, y);

				CGAL::Bounded_side PoincareSide = Poincare.bounded_side(t);
				CGAL::Bounded_side CircASide    = EuclidCircA.bounded_side(t);
				CGAL::Bounded_side CircBbSide   = EuclidCircBb.bounded_side(t);

				// First off, the point needs to be inside the Poincare disk. if not, there's no hope.
				if ( PoincareSide == CGAL::ON_BOUNDED_SIDE ) {
					
					// Inside the Poincare disk, but still outside the fundamental domain
					if ( CircASide  == CGAL::ON_BOUNDED_SIDE || 
							 CircBbSide == CGAL::ON_BOUNDED_SIDE   ) {
						return CGAL::ON_UNBOUNDED_SIDE;
					}

					// Inside the Poincare disk and inside the fundamental domain
					if ( CircASide  == CGAL::ON_UNBOUNDED_SIDE && 
							 CircBbSide == CGAL::ON_UNBOUNDED_SIDE ) {
						return CGAL::ON_BOUNDED_SIDE;
					} 

					// This is boundary, but we only consider the upper half. The lower half means outside.
					if (on_open_side) {
						return CGAL::ON_UNBOUNDED_SIDE;
					} else {
						return CGAL::ON_BOUNDED_SIDE;
					}

				} else {
					return CGAL::ON_UNBOUNDED_SIDE;
				}

			}

		};


		Side_of_fundamental_octagon
		side_of_fundamental_octagon_object() const {
			return Side_of_fundamental_octagon();
		}



}; // class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H








