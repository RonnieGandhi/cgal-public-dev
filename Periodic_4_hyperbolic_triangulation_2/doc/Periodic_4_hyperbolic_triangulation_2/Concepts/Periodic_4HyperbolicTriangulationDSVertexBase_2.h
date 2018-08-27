// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalModifBegin
\cgalRefines TriangulationDSVertexBase_2

A refinement of the concept `TriangulationDSVertexBase_2` that adds an interface for hyperbolic translations.
\cgalModifEnd

For periodic hyperbolic triangulations, the vertex base class needs to temporarily store a hyperbolic 
translation during the insertion process. 
\cgalModifBegin
A boolean flag indicates whether the face stores a translation or not. The value of the flag is automatically
set when storing or removing a translation.
\cgalModifEnd

\cgalHasModel `CGAL::Periodic_4_hyperbolic_triangulation_ds_vertex_base_2`

\sa `TriangulationDataStructure_2`
\sa `Periodic_4HyperbolicTriangulationDSFaceBase_2`

*/


class Periodic_4HyperbolicTriangulationDSVertexBase_2 {
public:

  /// \name Types
  /// @{
	
	/*!
  \cgalModifBegin
    (note: was inheriting from `Vb`, which is not defined here)
  \cgalModifEnd
  */
  typedef typename TriangulationDSVertexBase_2::Face_handle 		    
                                                      Face_handle;
  
  /*!
  \cgalModifBegin
    Must be the same as the point type `Periodic_4HyperbolicDelaunayTriangulationTraits_2::Point_2`
    defined by the geometric traits of the triangulation.
  \cgalModifEnd
  */
	typedef unspecified_type	    		                  Point;
	
  /*!
  \cgalModifBegin
    Must be the same as the translation type `Periodic_4HyperbolicDelaunayTriangulationTraits_2::Hyperbolic_translation`
    defined by the geometric traits of the triangulation.
  \cgalModifEnd
  */
  typedef unspecified_type                    			  Hyperbolic_translation;
  /// @}

  /// \name Creation
  /// @{
  /*!
    Default constructor.
  */
	Periodic_4HyperbolicTriangulationDSVertexBase_2();

  /*!
    Construct a vertex that stores the point `p`.
  */
	Periodic_4HyperbolicTriangulationDSVertexBase_2(const Point & p);

  /*!
    Constructs a vertex that stores the point `p` and is incident to the face `fh`.
  */
  Periodic_4HyperbolicTriangulationDSVertexBase_2(const Point & p, Face_handle fh);

  /*!
    Constructs a vertex that is incident to the face `fh`.
  */
	Periodic_4HyperbolicTriangulationDSVertexBase_2(const Face_handle& fh);
  /// @}


  /// \name Hyperbolic translations API
  /// @{
  /*!
    Stores the translation `tr` in the vertex, and sets a flag indicating that the vertex stores a translation.
  */
	void set_translation(const Hyperbolic_translation& tr);

  /*!
    Returns the translation stored in the vertex.
  */
	Hyperbolic_translation translation();

  /*!
    Removes the translation stored in the vertex, and resets the translation flag.
  */
	void clear_translation();

  /*!
    Returns the value of a flag, indicating whether the vertex stores a translation or not.
  */
	bool get_translation_flag();
  ///@}

};
