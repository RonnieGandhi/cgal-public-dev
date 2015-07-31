namespace CGAL {
/*!
\ingroup PkgVisibility_2Classes

\brief This class is a model of the concept `Visibility_2` can answer visibility queries within
a polygon that may have holes.

\details The class template comprises the implementation of the algorithm of Takao Asano and Tetsuo Asano \cite aaghi-vpsesp-85. The algorithm, as the name of the class template suggests, requires preprocessing. The preprocessing takes \f$ O(n^2)\f$ time and \f$ O(n^2)\f$ space, which reduces the query time to \f$O(n)\f$.

The main preprocessing step is the dualization of all vertices of the input arrangment into an arrangment of lines. 
Computing this arrangement takes \f$ O(n^2)\f$ time and \f$ O(n^2)\f$ space. 
Using this arrangment it is possible to retrieve the angular order of all vertices around 
a query point, which is one of the essential steps to achive linear query time. For more details we refer to \cite aaghi-vpsesp-85. 



\tparam Arrangement_2 is the type of input polygonal environment and output visibility polygon.

\tparam RegularizationCategory indicates whether the output should be regularized. It can be
specified by one of the following: ::Tag_true or ::Tag_false, where ::Tag_false is the default value.


\cgalModels `Visibility_2` 

\sa `CGAL::Simple_polygon_visibility_2<Arrangement_2, RegularizationCategory>`
\sa `CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationCategory>`
\sa `CGAL::Triangular_expansion_visibility_2<Arrangement_2, RegularizationCategory>`


*/
template <typename Arrangement_2, typename RegularizationCategory = Tag_false>
class Preprocessed_rotational_sweep_visibility_2 {
public:

/// \name Types 
/// @{

 /*!
  The type of the input arrangement.
  */
   typedef Arrangement_2 Arrangement_2;

 /*!
  The type of the output arrangement.
  */
   typedef Arrangement_2 Visibility_arrangement_2;

 /*! 
   The 2D point type used for the queries.  
 */ 
  typedef Arrangement_2::Point_2 Point_2; 

  /*!
   Face_const_handle type of the input arrangement.
   */
  typedef Arrangement_2::Face_const_handle Face_const_handle;

  /*!
   Halfedge_const_handle type of the input arrangement.
   */
  typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;



/// @}


/// \name Tags 
/// @{
  /*! 
    identifies whether the regularized visibility area is computed. 
  */
  typedef RegularizationCategory Regularization_category;
  
  /*! 
    identifies that the class supports general polygons (i.e.\ with holes). 
  */
  typedef ::Tag_true Supports_general_polygon_category; 

  /*! 
    identifies that the class supports general simple polygons. 
  */
  typedef ::Tag_true Supports_simple_polygon_category; 
/// @}



/// \name Constructors 
/// @{

/*!
Default constructor creates an empty `Preprocessed_rotational_sweep_visibility_2` object that is not
attached to any arrangement yet.
*/
Preprocessed_rotational_sweep_visibility_2();

/*! 
Constructs a `Preprocessed_rotational_sweep_visibility_2` object that is attached to `arr`.
*/ 
Preprocessed_rotational_sweep_visibility_2(const Arrangement_2& arr);

/// @}

/// \name functions 
/// @{

/*!
Returns whether an arrangement is attached to the visibility object
*/
  bool is_attached() const;

/*!
Attaches the given arrangement to the visibility object and applies preprocessing.
In case the object is already attached to another arrangement, 
the visibility object gets detached before being attached to `arr`.
*/
  void attach(const Arrangement_2& arr);

/*!
Detaches the arrangement from the visibility object it is currently attached to
*/
  void detach();

/*!
Access to the attached arrangement
*/
  const Arrangement_2& arr() const;

/*! 
Computes the visibility region of `q` in the
face `f` of the arrangement that is attached to the visibility object. 
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point
\param f is the face of the arrangement in which the visibility region is computed
\param out_arr is the output arrangement 
\pre `f` is a face of `arr()` and represents a valid polygon. 
\pre `q` is in the interior of the given face `f`
\return a handle to the face in `out_arr` that represents the visibility region
*/ 
  typename Visibility_arrangement_2::Face_handle compute_visibility(const Point_2& q, const Face_const_handle f, Visibility_arrangement_2& out_arr) const;

/*!
Computes the visibility region of `q` that is on `e`. If `q` is an interior point of `e`, the computed visibility region is restricted to the halfplane indicated by `e`. If `q` is an endpoint of `e`, the visibility region is restricted by `e` and its next.
The visibility region of `q` will be stored in `out_arr`.
\param q is the query point
\param e the halfedge on which `q` is located
\param out_arr is the output arrangement
\pre `e` is a halfedge of `arr()`
\pre `q` is on `e`
\pre `q` equals to `e->target()->point()` if `q` is an endpoint of `e`
\return a handle to the face in `out_arr` that represents the visibility region
*/
  typename Visibility_arrangement_2::Face_handle compute_visibility(const Point_2& q, const Halfedge_const_handle e, Visibility_arrangement_2& out_arr) const;

/// @}

}; /* end Visibility_2 */
}  /* namespace CGAL */
