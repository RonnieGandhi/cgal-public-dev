#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_H

#include <CGAL/Buildings/Utils/Level_of_detail_building_splitter.h>
#include <CGAL/Buildings/Utils/Level_of_detail_building_outliner.h>
#include <CGAL/Buildings/Utils/Level_of_detail_diagonalize_traits.h>

#include <CGAL/Buildings/Roofs/Roof/Level_of_detail_building_roof_fitter.h>
#include <CGAL/Buildings/Roofs/Roof/Level_of_detail_building_roof_estimator.h>
#include <CGAL/Buildings/Roofs/Roof/Level_of_detail_building_roof_face_validator.h>

#include <CGAL/Buildings/Roofs/Level_of_detail_building_roofs_creator.h>

#include <CGAL/Buildings/Roofs/Estimation/Level_of_detail_building_roof_estimator_box_strategy.h>
#include <CGAL/Buildings/Roofs/Estimation/Level_of_detail_building_roof_estimator_hull_strategy.h>
#include <CGAL/Buildings/Roofs/Estimation/Level_of_detail_building_roof_estimator_alpha_strategy.h>

#include <CGAL/Buildings/Envelope/Level_of_detail_building_envelope_input.h>
#include <CGAL/Buildings/Envelope/Level_of_detail_building_envelope_creator.h>
#include <CGAL/Buildings/Envelope/Level_of_detail_building_partition_input.h>
#include <CGAL/Buildings/Envelope/Level_of_detail_building_partition_creator.h>

#include <CGAL/Buildings/Envelope/Associaters/Level_of_detail_building_envelope_plane_associater.h>
#include <CGAL/Buildings/Envelope/Associaters/Level_of_detail_building_partition_naive_plane_associater.h>
#include <CGAL/Buildings/Envelope/Associaters/Level_of_detail_building_partition_vote_based_plane_associater.h>

#include <CGAL/Buildings/Kinetic/Level_of_detail_building_kinetic_partition_input_creator.h>
#include <CGAL/Buildings/Kinetic/Level_of_detail_building_kinetic_partition_output_creator.h>

#include <CGAL/Buildings/Cdt/Level_of_detail_building_cdt.h>
#include <CGAL/Buildings/Visibility/Level_of_detail_buildings_visibility_3.h>

#include <CGAL/Buildings/Utils/Level_of_detail_local_mesh_builder.h>
#include <CGAL/Buildings/Utils/Level_of_detail_planar_region_merger.h>
#include <CGAL/Buildings/Utils/Level_of_detail_buildings_facets_cleaner_3.h>

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_H