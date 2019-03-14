#include <boost/variant.hpp>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_3D_off_io.h>
#include <CGAL/config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/iterator.h>
 // For CGAL < 4.11
#if CGAL_VERSION_NR < 1041100000
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#endif  // CGAL_VERSION_NR < 1041100000
#include <fstream>
#include <cmath>
#include <string>
#include <tuple>
#include <map>
#include <utility>
#include <vector>
#include <cstdlib>
#include <utilities/Alpha_complex/alpha_complex_3d_helper.h>
// Alpha_shape_3 templates type definitions
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
// For CGAL < 4.11
#if CGAL_VERSION_NR < 1041100000
using Gt = CGAL::Regular_triangulation_euclidean_traits_3<Kernel>;
using Vb = CGAL::Alpha_shape_vertex_base_3<Gt>;
using Fb = CGAL::Alpha_shape_cell_base_3<Gt>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Fb>;
using Triangulation_3 = CGAL::Regular_triangulation_3<Gt, Tds>;
// From file type definition
using Point_3 = Gt::Bare_point;
using Weighted_point_3 = Gt::Weighted_point;
// For CGAL >= 4.11
#else   // CGAL_VERSION_NR < 1041100000
using Rvb = CGAL::Regular_triangulation_vertex_base_3<Kernel>;
using Vb = CGAL::Alpha_shape_vertex_base_3<Kernel, Rvb>;
using Rcb = CGAL::Regular_triangulation_cell_base_3<Kernel>;
using Cb = CGAL::Alpha_shape_cell_base_3<Kernel, Rcb>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Triangulation_3 = CGAL::Regular_triangulation_3<Kernel, Tds>;
// From file type definition
using Point_3 = Triangulation_3::Bare_point;
using Weighted_point_3 = Triangulation_3::Weighted_point;
#endif  // CGAL_VERSION_NR < 1041100000
using Alpha_shape_3 = CGAL::Alpha_shape_3<Triangulation_3>;
// filtration with alpha values needed type definition
using Alpha_value_type = Alpha_shape_3::FT;
using Object = CGAL::Object;
using Dispatch =
CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<Object, Alpha_value_type>,
	CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<Object> >,
	std::back_insert_iterator<std::vector<Alpha_value_type> > > >;
using Cell_handle = Alpha_shape_3::Cell_handle;
using Facet = Alpha_shape_3::Facet;
using Edge_3 = Alpha_shape_3::Edge;
using Vertex_handle = Alpha_shape_3::Vertex_handle;
using Vertex_list = std::vector<Alpha_shape_3::Vertex_handle>;
// gudhi type definition
using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = ST::Filtration_value;
using Simplex_tree_vertex = ST::Vertex_handle;
using Alpha_shape_simplex_tree_map = std::map<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex>;
using Simplex_tree_vector_vertex = std::vector<Simplex_tree_vertex>;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<ST, Gudhi::persistent_cohomology::Field_Zp>;

struct cmp_intervals_by_dim_then_length {
	explicit cmp_intervals_by_dim_then_length(ST * sc)
		: sc_(sc) { }
	template<typename Persistent_interval>
	bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
		if (sc_->dimension(get < 0 >(p1)) == sc_->dimension(get < 0 >(p2)))
			return (sc_->filtration(get < 1 >(p1)) - sc_->filtration(get < 0 >(p1))
> sc_->filtration(get < 1 >(p2)) - sc_->filtration(get < 0 >(p2)));
		else
			return (sc_->dimension(get < 0 >(p1)) > sc_->dimension(get < 0 >(p2)));
	}
	ST* sc_;
};