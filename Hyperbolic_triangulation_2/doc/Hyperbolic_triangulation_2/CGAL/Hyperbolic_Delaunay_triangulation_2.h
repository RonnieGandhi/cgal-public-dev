// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.


namespace CGAL {
 
/*!
\ingroup PkgHyperbolicTriangulation2MainClasses

The class `Hyperbolic_Delaunay_triangulation_2` is the main class of the 2D Hyperbolic Triangulations package.
It is designed to represent Delaunay triangulations of sets of points in the hyperbolic plane. 
The hyperbolic plane is represented in the Poincaré disk model.

\tparam Gt  must be a model of `HyperbolicDelaunayTriangulationTraits_2`.
\tparam Tds must be a model of `TriangulationDataStructure_2`. By default, this parameter is instantiated with 
`Triangulation_data_structure_2< Triangulation_vertex_base_2<Gt>, Hyperbolic_triangulation_face_base_2<Gt> >`

\sa HyperbolicDelaunayTriangulationTraits_2
\sa TriangulationDataStructure_2
\sa Delaunay_triangulation_2

*/

template < class Gt, class Tds >
class Hyperbolic_Delaunay_triangulation_2: private Delaunay_triangulation_2<Gt,Tds> {

public:

  /// \name Types
  /// @{
    typedef Gt                                          Geom_traits;
    typedef Tds                                         Triangulation_data_structure;
  /// \name 
  

  /// \name
  /// @{
    /*!
      Size type (integral unsigned).
    */
    typedef typename Triangulation_data_structure::size_type        size_type;
    typedef typename Geom_traits::Point_2                           Point;
    typedef typename Geom_traits::Triangle_2                        Triangle;
  /// @}


  /// \name
  /// The following types are defined to give access to the elements of the triangulation:
  /// @{
    typedef typename Triangulation_data_structure::Vertex_handle    Vertex_handle;
    typedef typename Triangulation_data_structure::Face_handle      Face_handle;
    typedef typename Triangulation_data_structure::Edge             Edge;
  /// @}


  /// \name 
  /// The following types are defined for use in the construction of the Voronoi diagram:
  /// @{
    typedef typename Geom_traits::Voronoi_point_2       Voronoi_point;
    typedef typename Geom_traits::Hyperbolic_segment_2  Hyperbolic_segment;
  /// @}

    
  /// \name 
  /// \cgalModifBegin
  /// The following iterator and circulator types are defined to give access over hyperbolic faces and edges:
  /// \cgalModifEnd
  /// @{
    typedef unspecified_type                            All_faces_iterator;
    typedef unspecified_type                            All_edges_iterator;
    typedef unspecified_type                            All_vertices_iterator;
    typedef unspecified_type                            Vertex_circulator;
  /// @}

  /// \name
  /// \cgalModifBegin
  /// The enumeration `Locate_type` is defined to specify which case occurs when locating a 
  /// point in the triangulation.
  /// \cgalModifEnd
  /// @{
    enum Locate_type { 
      VERTEX = 0, 
      EDGE, 
      FACE, 
      OUTSIDE_CONVEX_HULL, 
      OUTSIDE_AFFINE_HULL 
    };
  /// @}  

  /// \name Creation
  /// @{    
    /*!
      %Default constructor
    */
    Hyperbolic_Delaunay_triangulation_2(const Geom_traits& gt = Geom_traits());
    
    /*!
      Copy constructor
    */
    Hyperbolic_Delaunay_triangulation_2(const Hyperbolic_Delaunay_triangulation_2<Gt,Tds> &tr);

    /*!
      Equivalent to creating an empty triangulation and calling `insert(first, last)`.
    */
    template<class InputIterator>
    Hyperbolic_Delaunay_triangulation_2(InputIterator first, InputIterator last, const Geom_traits& gt = Geom_traits());
  /// @}


  /// \name Assignment
  /// @{

    /*!
      The triangulation `tr` is duplicated, and modifying the copy after the duplication does not modify the original.
      \todo implement! 
    */
    Hyperbolic_Delaunay_triangulation_2& operator=(const Hyperbolic_Delaunay_triangulation_2& tr);
    
    /*!
      The triangulation is swapped with `tr`.
      \todo implement!
    */
    void swap(Hyperbolic_Delaunay_triangulation_2& tr);
    
    /*!
      Deletes all finite vertices and faces of the triangulation.
      \todo implement!
    */
    void clear();

    /*!
      Equality operator.
      \todo implement!
    */
    bool operator==(const Hyperbolic_Delaunay_triangulation_2<Gt,Tds>& t1, 
                    const Hyperbolic_Delaunay_triangulation_2<Gt,Tds>& t2);
    
    /*!
      Inequality operator.
      \todo implement!
    */
    bool operator!=(const Hyperbolic_Delaunay_triangulation_2<Gt,Tds>& t1, 
                    const Hyperbolic_Delaunay_triangulation_2<Gt,Tds>& t2);
  /// @}


  /// \name Access functions
  /// @{
    /*!
      Returns a const reference to the geometric traits object.
      \todo implement!
    */
    const Geom_traits&    geom_traits() const;

    /*!
      Returns a const reference to the triangulation data structure.
      \todo implement!
    */
    const Triangulation_data_structure& tds() const; 

    /*!
      Returns a reference to the triangulation data structure.
      \todo implement!
    */
    Triangulation_data_structure& tds();

    /*!
      Returns the dimension of the affine hull.
    */
    int dimension() const;

    /*!
      Returns the number of vertices.
    */
    size_type             number_of_vertices()          const;
    
    /*!
      Returns the number of hyperbolic edges.
    */
    size_type             number_of_hyperbolic_edges()  const;
    
    /*!
      Returns the number of hyperbolic faces.
    */
    size_type             number_of_hyperbolic_faces()  const;
  /// @}


  /// \name Geometric access functions
  /// @{
    /*
      Returns the triangle formed by the three vertices of face `f`.
    */
    Triangle triangle(Face_handle f) const;

    /*!
      \cgalModifBegin
      Returns the hyperbolic segment formed by the vertices of the edge `(f, i)`.
      \cgalModifEnd
    */
    Hyperbolic_segment segment(Face_handle f, int i) const;

    /*!
      \cgalModifBegin
      Returns the hyperbolic segment formed by the vertices of edge `e`.
      \cgalModifEnd
    */
    Hyperbolic_segment segment (const Edge& e) const;
  ///@}


  /// \name Insertion
  /// @{
    /*!
      Inserts the point `p` in the triangulation.
      If the point `p` coincides with a existing vertex, then the vertex is returned 
      and the triangulation is not modified. The optional parameter `start` is used
      to initialize the location of `p`.
    */
    Vertex_handle insert(const Point  &p, Face_handle start = Face_handle() );
    
    /*!
      Inserts the point p at the location given by `(lt,loc,li)`. 
      The handle to the new vertex is returned.

      \sa locate()
    */
    Vertex_handle insert(const Point& p, typename Locate_type lt, Face_handle loc, int li );
      

    /*!
      Inserts the points in the range [first,last) into the triangulation.
      Returns the number of inserted points. Note that this function is not 
      guaranteed to insert the points following the order of `InputIterator`, 
      as `spatial_sort()` is used to improve efficiency.

      \tparam InputIterator  must be an input iterator with the value type `Point_2`.
    */
    template < class InputIterator >
    std::ptrdiff_t insert(InputIterator first, InputIterator last);
  /// @}

  /// \name Removal
  /// @{
    /*!
      Removes the vertex `v` from the triangulation.
      \pre `v` is a vertex of the triangulation.
      \todo implement this function!
    */
    void remove(Vertex_handle v);
  /// @}


  /// \name Point Location
  /// @{
    /*!
      Locates the point `query` in the triangulation.

      If the point `query` lies inside the hyperbolic convex hull of the points of the triangulation,
      then the hyperbolic face that contains the query in its interior is returned. 

      If `query` lies on a vertex or on an edge, then one of the faces having `query` on its boundary is returned.
      \todo verify case when `query` lies on dangling edge!
    
      If `query` lies outside of the convex hull of the points of the triangulation, then a default-constructed
      `Face_handle()` is returned. 

      The optional argument `hint` is used as a starting place for the search.
    */
    Face_handle locate(const Point& query, const Face_handle hint = Face_handle()) const;


    /*!
      Same as above. 

      The variable `lt` contains information about the element in which `query` has been located. See 
      the enumeration `Triangulation::Locate_type` for details.

      If `lt` is `Triangulation::Locate_type::VERTEX`, then the variable `li` contains the index of the 
      vertex in the returned face. If `lt` is `Triangulation::Locate_type::EDGE`, then `li` is the index 
      of the edge in the returned face. 
    */
    Face_handle locate(const Point& query, Locate_type& lt, int &li, Face_handle hint = Face_handle()) const;
  /// @}


  /// \name Queries
  /// @{

    /*!
      Computes the conflict zone induced by `p`.

      If the optional parameter `start` is given, then it must be a face in conflict with `p`. 
      Returns an iterator on the faces of the triangulation in conflict with `p`.
    */
    template<class OutputItFaces>
    OutputItFaces find_conflicts(const Point& p, 
                                OutputItFaces fit, 
                                Face_handle start = Face_handle()) const;
  ///@}
  

  /// \name Vertex, Face and Edge iterators and circulators
  /// @{
    
    All_vertices_iterator   all_vertices_begin()  const;
    All_vertices_iterator   all_vertices_end()    const;

    All_edges_iterator      all_edges_begin() const;
    All_edges_iterator      all_edges_end()   const;
    
    All_faces_iterator      all_faces_begin() const;
    All_faces_iterator      all_faces_end()   const;  

    Vertex_circulator       adjacent_vertices(Vertex_handle v) const;
  /// @}


  
  /// \name Voronoi Diagram 
  /// Note that the user should use a kernel with exact constructions in order to guarantee 
  /// the computation of the Voronoi diagram (as opposed to computing the triangulation only, 
  /// which requires only exact predicates).
  /// @{
    /*!
      Returns the hyperbolic center of the circumdisk of `f`.
      \pre `f` is hyperbolic
    */
    Voronoi_point dual(Face_handle f) const;

    /*!
      Returns the hyperbolic segment that is dual to `e`.
    */
    Hyperbolic_segment dual(const Edge& e) const;

    /*!
      Returns the hyperbolic segment that is dual to the edge `(f,i)`.
      \pre `f` is hyperbolic
    */
    Hyperbolic_segment dual(Face_handle f, int i) const;
  /// @}

};
  
} //namespace CGAL

