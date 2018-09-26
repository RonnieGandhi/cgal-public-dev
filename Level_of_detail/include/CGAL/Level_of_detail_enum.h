#ifndef CGAL_LEVEL_OF_DETAIL_ENUM_H
#define CGAL_LEVEL_OF_DETAIL_ENUM_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#define PN "\r\n"
#else 
#define PSR "/" 
#define PN "\n"
#endif

// STL includes.
#include <map>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CGAL {

	namespace LOD {

		enum class Structured_label { 
			LINEAR, 
			CORNER, 
			CLUTTER 
		};

		// Visibility labels.
		enum class Visibility_label { 
			IN,  	// inside a triangle
			OUT, 	// outside a triangle
			UNKNOWN // 50/50 inside outside
		};

		// Visibility methods.
		enum class Visibility_method {
			POINT_BASED_CLASSIFICATION,    // this is a repeatable method
			FACE_BASED_BARYCENTRIC,        // does not work
			FACE_BASED_NATURAL_NEIGHBOURS, // every time can give different result since it is randomized
			FACE_BASED_AFFINE,			   // use affine coordinates instead of natural neighbours
			FACE_BASED_COUNT			   // hard to find the proper circle radius
		};

		// Visibility approaches.
		enum class Visibility_approach {
			POINT_BASED, // here we go over all input points
			FACE_BASED   // here we go over all input faces
		};

		// Visibility samplers.
		enum class Visibility_sampler {
			RANDOM_UNIFORM_0,    // a bit faster than UNIFORM_1 but the result is similar, randomized
			RANDOM_UNIFORM_1,    // see above
			UNIFORM_SUBDIVISION, // determenistic sampler based on midpoint subdivision
			BARYCENTRE			 // barycentre of the given triangle
		};

		// Type of the building's boudnary.
		enum class Building_boundary_type {
			ORIENTED,  // clockwise-oriented boundary
			UNORIENTED // not oriented boundary, which is the set of segments
		};

		// Face info class.
		template<typename FT>
		class My_face_info {

		public:
			FT in = FT(1) / FT(2); 				  			 // visibility label (in - inside) or (out - outside); or alternatively label A - everything that is bigger than 0.5, label B < 0.5
			CGAL::Color in_color = CGAL::Color(255, 205, 0); // visibility color (in < 1/2 - red), (in = 1/2 - yellow), (in > 1/2 - green)
			
			int bu = -1; 		   						 	   // building's index - (0, 1, 2 etc.) where (-1 means not a building)
			CGAL::Color bu_color = CGAL::Color(169, 169, 169); // building's color - random color is used per building

			bool is_valid = true; // check if this is a valid face

			CGAL::Color color = CGAL::Color(255, 255, 255); // default face color
			
			int index 	   = -1;
			int roof_index = -1;

			CGAL::Color tag_color = CGAL::Color(255, 205, 0);
		};

		// Vertex info class.
		template<typename Label, typename FT>
		class My_vertex_info {

		public:
			Label label = Label::CLUTTER;
			CGAL::Color color = CGAL::Color(0, 0, 0);
			FT height = -FT(100000000000000);
			int index = -1;
		};

		// Building structure.
		template<class Kernel, class InputCDT>
		struct Building {

		public:
			using FT  		 = typename Kernel::FT;
			using Point_3    = typename Kernel::Point_3;
			using Triangle_3 = typename Kernel::Triangle_3;
			using Segment_2  = typename Kernel::Segment_2;
			using Segment_3  = typename Kernel::Segment_3;
			using Plane_3    = typename Kernel::Plane_3;

			using CDT = InputCDT;

			using Face_handle   = typename CDT::Face_handle;
			using Vertex_handle = typename CDT::Vertex_handle;

			FT height 		  = FT(0); 				  // height of the building
			CGAL::Color color = CGAL::Color(0, 0, 0); // color of the building

			FT roofs_min_height = FT(0); // min height among all roof points
			FT roofs_max_height = FT(0); // max height among all roof points

			FT max_height = -FT(1);

			struct Roof {
				
				using Roof_boundary     = std::vector<Point_3>;
				using Associated_planes = std::vector<size_t>;

				Roof_boundary     boundary;
				Associated_planes associated_planes;

				bool is_plane_index = false;

				bool is_valid = true; // debugging info
				Roof_boundary tmp; 	  // not needed in the final version!

				int index = -1;
			};

			struct Data {
                size_t index;

				bool is_vertical = false;
                CGAL::Color color;
            };

			struct Wall {
				using Wall_boundary = std::vector<Point_3>;
				Wall_boundary boundary;
			};

			using Walls = std::vector<Wall>;
			Walls walls;

            using Data_triangle  = std::pair<Triangle_3, Data>;
            using Data_triangles = std::vector<Data_triangle>;

			struct Partition_element {

				Segment_3 segment;
				bool to_be_used = false;
			};

			using Partition_input    = std::vector<Partition_element>;
			using Partition_segments = std::vector<Segment_2>;

			using Index   	 = int;
			using Indices 	 = std::vector<Index>;
			using Shapes  	 = std::vector<Indices>;
			using Roofs      = std::vector<Roof>;
			using Boundary 	 = std::vector<Vertex_handle>;
			using Boundaries = std::vector<Boundary>;
			using Planes 	 = std::vector<Plane_3>;

			using Floor_faces = std::vector<Face_handle>;

			Boundaries 		   				  								  boundaries; // boundary vertices of the building ordered counterclockwise (may store multiple boundaries)
			std::vector< std::map<Vertex_handle, std::vector<Face_handle> > > wedges;     // all faces adjacent to each boundary vertex above - must be unique face handles
			Floor_faces   						   				  			  faces;	  // all faces that belong to this building

			bool is_oriented = true;    		// flag to check if the computed boundary is oriented or not, see Building_boundary_type above
			std::unordered_set<int> neighbours; // indices of all neighbouring buildings of the given building

			Indices 	       interior_indices;   // indices of all input points that lie inside this building
			Shapes  	       shapes; 		       // detected shapes by region growing
			Roofs   	       roofs;			   // roofs = bounding boxes of points projected on the respected planes found in shapes above	
			Data_triangles     envelope_input;     // input for the 3D envelope
			Partition_input    partition_input;    // input for the roof partitioning
			Planes 			   planes;			   // all roof planes associated with this building
			Partition_segments partition_segments; // 2D segments used to compute roof partition

			CDT cdt; 			  // cdt used for roofs
			bool is_valid = true; // flag to check if we should output this building or not, if it is a valid building or not
			int index = -1;

            using Contribution  		 	= int;
            using Contributions 		 	= std::vector<Contribution>;
            using Face_contributions 	 	= std::map<Face_handle, Contributions>;
			using Contribution_pair      	= std::pair<Face_handle, Contributions>;
            using Comparator             	= std::function<bool(Contribution_pair, Contribution_pair)>;
            using Sorted_face_contributions = std::vector<Contribution_pair>;

			Face_contributions 	   	  face_contributions;
			Sorted_face_contributions sorted_face_contributions;

			FT 	   current_percentage;
			size_t total_contributions_size;

			struct Graphcut_facet {
			
			public:
				using Data 		= std::pair<int, int>;
				using Data_pair = std::pair<Data, Data>;

				Data 	  data;
				Data_pair neighbours;

				FT weight  = -FT(1);
				FT quality = -FT(1);

				bool is_valid = true;
			};

			using Graphcut_facets = std::vector<Graphcut_facet>;

			Graphcut_facets graphcut_facets;

			struct Polyhedron {
			
			public:
				using Vertex   = Point_3;
				using Vertices = std::vector<Vertex>;
				
				struct Facet {
					using Indices = std::vector<int>;
					
					Indices indices;
					bool is_valid = true;
				};

				using Facets = std::vector<Facet>;

				Vertices vertices;
				Facets 	 facets;

				bool is_valid = true;

				FT in  = -FT(1);
				FT out = -FT(1);

				FT weight = -FT(1);
			};

			using Polyhedrons = std::vector<Polyhedron>;
			
			Polyhedrons polyhedrons;
			Polyhedrons polyhedron_facets;

			using JP_point_3  = typename CGAL::Exact_predicates_exact_constructions_kernel::Point_3;
            using JP_polygon  = std::vector<JP_point_3>;
            using JP_polygons = std::vector<JP_polygon>;

			JP_polygons jp_polygons;

			bool is_clean = false;

			using Clean_facet  = std::vector<Point_3>;
			using Clean_facets = std::vector<Clean_facet>;
			
			Clean_facets clean_facets;

			using Region_facet   = std::vector<Point_3>;
            using Region_facets  = std::vector<Region_facet>;
			using Output_regions = std::vector<Region_facets>;

			Output_regions output_regions;

			void clear_interior_indices() {
				interior_indices.clear();
			}

			void clear_shapes() {
				shapes.clear();
			}

			void clear_roofs() {
				roofs.clear();
			}

			void clear_envelope_input() {
				envelope_input.clear();
			}

			void clear_partition_input() {
				partition_input.clear();
			}

			void clear_planes() {
				planes.clear();
			}
			
			void clear_partition_segments() {
				partition_segments.clear();
			}
		};

		// Type of the roof fitter.
		enum class Roof_fitter_type {
			MIN, // fit data to the minimum height
			AVG, // fit data to the average height
			MAX  // fit data to the maximum height
		};

		// Main test data.
		enum class Main_test_data_type {
			BASIC,           // basic data set from the Loader_stub class.
			COMPLEX,         // sketch up generated simple data set with square buildings
			PARIS,           // half of the paris real data set
			P10,             // p10 data set from the original LOD paper of Yannick
			PARIS_FULL,      // full paris data set
			PARIS_ETH,       // half paris data set classified with ETH random forest
			PARIS_FULL_ETH,  // full paris data set classified with ETH random forest
			RESIDENT_TILE_1, // three different residential tiles
			RESIDENT_TILE_2,
			RESIDENT_TILE_3,
			PARIS_TILE_1, 	 // two different Paris tiles
			PARIS_TILE_2,
			PARIS_BIG        // 1 km - 9 tiles stitched together
		};

		// Nearest neighbour search.
		enum class Neighbour_search_type {
			KNN,    // use k nearest neighbours
			CIRCLE, // use all points from the circle of the given radius
			SQUARE  // use all points from the square of the given radius
		};

		// Fitter type used in thinning.
		enum class Thinning_fitter_type { 
			LINE // fit to the line
		};

		// New point type used in the grid simplify algorithm.
		enum class Grid_new_point_type { 
			CENTROID,   // centroid of the cell
			BARYCENTRE, // barycentre of the given set of samples
			CLOSEST     // point closest to the barycentre above
		};

		// Thinning type.
		enum class Thinning_type {
			NAIVE,	// naive thinning where we project all points onto a line
			COMPLEX // a complex version, where we perform many different optimization steps and preserve features
		};

		// Thinning scale type.
		enum class Thinning_scale_type {
			FIXED,			// fixed manually set scale
			ADAPTIVE_FIXED, // automatically chosen scale
			PROGRESSIVE		// scale changes progressively
		};

		// Structuring corner algorithm.
		enum class Structuring_corner_algorithm {
			NO_CORNERS,         // we do not insert any corners
			GRAPH_BASED,	    // in this algorithm, we build an adjacency graph and insert corners based on this graph
			INTERSECTION_BASED, // in this algorithm, we intersect all segments and insert the best intersections
			NO_T_CORNERS		// we do not insert T-like corners
		};

		// Structuring adjacency threshold method.
		enum class Structuring_adjacency_threshold_method {
			LOCAL, // internal local epsilon is chosen
			GLOBAL // user-defined value is chosen
		};

		// Method to estimate normals in 2D region growing.
		enum class Region_growing_normal_estimation {
			PROJECTED, // project exact normals if they exist
			LOCAL      // estimate normals using PCA
		};

		// Quality data type.
		enum class Quality_data_type { 	
			DST_MIN,   // distortion types
			DST_AVG, 
			DST_MAX,
			DST_ROOFS,
			DST_WALLS,
			CMP_ROOFS, // complexity metric
			CMP_WALLS,
			CMP,
			COV_ROOFS, // coverage metric
			COV_WALLS,
			COV,       
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUM_H