// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgHyperbolicTriangulation2Concepts
\cgalConcept

The concept `HyperbolicDelaunayTriangulationTraits_2` describes the set of requirements 
to be fulfilled by any class used to instantiate the first template parameter of the class
`CGAL::Hyperbolic_Delaunay_triangulation_2<Traits, Tds>`. It defines the geometric objects 
(points, segments...) forming the triangulation together with geometric predicates and 
constructions on these objects.

\cgalHasModel CGAL::Hyperbolic_Delaunay_triangulation_traits_2
\cgalHasModel CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2
*/


class HyperbolicDelaunayTriangulationTraits_2 {


public:

  /// \name Types
  /// @{

	/*!
		Represents a point in the Poincaré disk or on the (unit) circle at infinity.
	*/
	typedef unspecified_type 			Point_2;

	/*!
		Represents the dual object of a triangle in the hyperbolic Delaunay triangulation. 
		The dual of a Delaunay triangle is the <i>hyperbolic</i> center of the circle circumscribing it.
	*/
	typedef unspecified_type 		    Voronoi_point_2;

	/*!
		Represents a hyperbolic segment defined by two points. 
		In the Poincaré disk model, a hyperbolic segment is supported either by the Euclidean circle 
		that passes through the two points and is perpendicular to the circle at infinity, or by the 
		Euclidean line that passes through the two points and the origin. Abusively, we allow one or 
		both endpoints of the segment to lie on the circle at infinity, so a hyperbolic segment can 
		actually represent a hyperbolic ray or a hyperbolic line.
	*/
	typedef unspecified_type	        Hyperbolic_segment_2;

	/*!
		Represents a triangle in the hyperbolic plane defined by three points. 
	*/
	typedef unspecified_type 			Triangle_2;
  /// @}

  /// \name Predicate Types
  /// @{
	/*!
		A predicate object. Must provide the function operator

		`Orientation operator()(Point_2 p, Point_2 q, Point_2 r),`
		
		which returns the orientation of the points `p, q`, and `r`.
	*/
	typedef unspecified_type			Orientation_2;

	/*!
		A predicate object. Must provide the function operator

		`Oriented_side operator()(Point_2 p, Point_2 q, Point_2 r, Point_2 t),`
		
		which returns the position of point `t` relative to the oriented circle
		defined by the points `p, q`, and `r`. 
	*/
	typedef unspecified_type			Side_of_oriented_circle_2;


	/*!
		\cgalModifBegin
		A predicate object. Must provide the function operator

		`Oriented_side operator()(Point_2 p, Point_2 q, Point_2 query),`

		which returns the position of point `query` relative to the oriented hyperbolic 
		segment with vertices `p` and `q`. 
		\cgalModifEnd
	*/
	typedef unspecified_type 			Side_of_oriented_hyperbolic_segment_2;



  /// @}

  /// \name Construction Types
  /// @{
	/*!
		A constructor object. 

		Must provide the function operator
		
		`Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q),`
		
		which constructs a hyperbolic segment from two points `p` and `q`.
		Note that `p` and `q` may also lie on the circle at infinity.
	*/
	typedef unspecified_type 			Construct_hyperbolic_segment_2;

	/*!
		A constructor object. Must provide the function operator

		`Voronoi_point_2 operator()(Point_2 p, Point_2 q, Point_2 r),`
		
		which constructs the hyperbolic circumcenter of the triangle with 
		vertices `p, q`, and `r`.
	*/
	typedef unspecified_type 			Construct_hyperbolic_circumcenter_2;

	/*!
		A constructor object. Must provide the function operator

		`Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q),`
		
		which constructs the hyperbolic bisector of two points `p` and `q` lying 
		in the Poincaré disk. The endpoints of the resulting hyperbolic segment 
		lie on the circle at infinity. It must also provide the function operator

		`Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q, Point_2 r),`

		where the points `p, q`, and `r` lie in the Poincaré disk. This overloaded 
		version constructs the hyperbolic bisector of the segment [p,q] limited by 
		the hyperbolic circumcenter of `p, q, r` on one side and the circle at 
		infinity on the other. Moreover, it must provide the function operator

		`Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q, Point_2 r, Point_2 s),`

		where the points `p, q, r`, and `s` lie in the Poincaré disk. This overloaded
		version constructs the hyperbolic bisector of the segment [p,q] limited by 
     	the hyperbolic circumcenter of `p, q, r` on one side, and the hyperbolic 
     	circumcenter of `p, s, q` on the other side.
	*/
	typedef unspecified_type 			Construct_hyperbolic_bisector_2;
  /// @}


  /// \name Operations
  /// The following functions give access to the predicate objects.
  /// @{  
	Orientation_2                          orientation_2_object();
	Side_of_oriented_circle_2              side_of_oriented_circle_2_object();  
	/*!
		\cgalModifBegin
	*/
	Side_of_oriented_hyperbolic_segment_2  side_of_oriented_hyperbolic_segment_2_object();
	/*!
		\cgalModifEnd
	*/  
  /// @}

  /// \name
  /// The following functions must be provided only if the methods of `Hyperbolic_Delaunay_triangulation_2`
  /// that return elements of the Voronoi diagram are instantiated:
  /// @{
	Construct_hyperbolic_segment_2       construct_hyperbolic_segment_2_object();
	Construct_hyperbolic_circumcenter_2  construct_hyperbolic_circumcenter_2_object();
	Construct_hyperbolic_bisector_2      construct_hyperbolic_bisector_2_object();
  /// @}
  

}; 

