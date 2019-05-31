#pragma once

//-- Boost Includes
#include <boost/variant.hpp>

//-- GUDHI Includes
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_3D_off_io.h>
//#include <gudhi/Persistence_graph.h>
//#include <gudhi/Persistence_intervals.h>
//#include <gudhi/Persistence_heat_maps.h>

//-- Local Includes
#include "cgal_helper.h"

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = ST::Filtration_value;
using Simplex_tree_vertex = ST::Vertex_handle;
using Alpha_shape_simplex_tree_map = std::map<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex>;
using Simplex_tree_vector_vertex = std::vector<Simplex_tree_vertex>;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<ST, Gudhi::persistent_cohomology::Field_Zp>;
//using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;
//using Persistence_heatmap = Gudhi::Persistence_representations::Persistence_heat_maps<Gudhi::Persistence_representations::constant_scaling_function>;
