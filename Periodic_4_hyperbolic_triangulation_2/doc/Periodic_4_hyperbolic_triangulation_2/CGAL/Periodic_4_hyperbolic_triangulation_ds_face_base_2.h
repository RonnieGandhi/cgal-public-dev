// Copyright (c) 1999-2018   INRIA Nancy - Grand Est (France).
// All rights reserved.


namespace CGAL {

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2VertexFaceClasses

The class `Periodic_4_hyperbolic_triangulation_ds_face_base_2` is the default model for the 
concept `Periodic_4HyperbolicTriangulationDSFaceBase_2`. It accepts two template parameters:

\tparam GT 	Geometric Traits type. This should be a model of the concept `Periodic_4HyperbolicDelaunayTriangulationTraits_2`.
			This template parameter has no default value.
\tparam FB 	Face Base type. Should be a model of the concept `TriangulationDSFaceBase_2`. 
		   	The default value for this template parameter is `Triangulation_ds_face_base_2<>`

\cgalModels Periodic_4HyperbolicTriangulationDSFaceBase_2

\sa `Periodic_4_hyperbolic_triangulation_ds_vertex_base_2` 
*/



template< typename GT, typename FB >
class Periodic_4_hyperbolic_triangulation_ds_face_base_2 : public FB {

};



}  // namespace CGAL

