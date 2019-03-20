#include <boost/version.hpp>
#include <boost/program_options.hpp>
#include <boost/variant.hpp>
#include <boost/algorithm/string.hpp>

#include "cifio.h"
#include "math.h"  // NOLINT(modernize-deprecated-headers)
#include "file.h"

#include <muParser.h>

#include "atom_data.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include "persistence_helper.h"

#include "niggli.h"

struct crystal_structure
{
	std::vector<std::string> atom_names;
	std::vector<Eigen::Vector3d> atom_vectors;
	struct
	{
		double a, b, c, alpha, beta, gamma, volume;
	} unit_cell;

	void clear() {
		atom_names.clear();
		atom_vectors.clear();

		unit_cell.a = unit_cell.b = unit_cell.c = unit_cell.alpha = unit_cell.beta = unit_cell.gamma = unit_cell.volume = 0;
	}
};

struct parser_parameters
{
	bool show_hydrogen;
	bool trim_unit_cell;
	bool reduce_unit_cell;
	bool check_duplicates;
	bool offset_by_min;
	unsigned int extension;
};

void parse_file(const std::string& filepath, CIFHandler& ch, crystal_structure& cs, parser_parameters& pp) {
	auto root = ch.copyCIFrom(cif::read_file(filepath))->blocks[0];

	cs.clear();

	cs.unit_cell.a = io::stod(*root.find_value("_cell_length_a"));
	cs.unit_cell.b = io::stod(*root.find_value("_cell_length_b"));
	cs.unit_cell.c = io::stod(*root.find_value("_cell_length_c"));
	cs.unit_cell.alpha = io::stod(*root.find_value("_cell_angle_alpha"));
	cs.unit_cell.beta  = io::stod(*root.find_value("_cell_angle_beta"));
	cs.unit_cell.gamma = io::stod(*root.find_value("_cell_angle_gamma"));

	if (pp.reduce_unit_cell) {
		std::vector<double> reduced_values = {
			cs.unit_cell.a,
			cs.unit_cell.b,
			cs.unit_cell.b,
			cs.unit_cell.alpha,
			cs.unit_cell.beta,
			cs.unit_cell.gamma
		};

		niggli::reduce(reduced_values, 0.0001);

		cs.unit_cell.a = reduced_values[0];
		cs.unit_cell.b = reduced_values[1];
		cs.unit_cell.c = reduced_values[2];
		cs.unit_cell.alpha = reduced_values[3];
		cs.unit_cell.beta = reduced_values[4];
		cs.unit_cell.gamma = reduced_values[5];
	}

	const auto atom_loop = root.find_loop("_atom_site_label").get_loop();

	const int type_pos = atom_loop->find_tag("_atom_site_type_symbol");
	const int xyz_pos  = atom_loop->find_tag("_atom_site_fract_x");
	

	const int val_size = atom_loop->values.size();
	const int tag_size = atom_loop->tags.size();

	const int atoms_size = val_size / tag_size;

	cs.atom_vectors.reserve(atoms_size);
	cs.atom_names.reserve(atoms_size);

	std::string sym_string = "_symmetry_equiv_pos_as_xyz";
	auto sym_tag = root.find_loop("_symmetry_equiv_pos_as_xyz");
	
	gemmi::cif::Loop* sym_loop = sym_tag.get_loop();
	int sym_pos = -1;
	
	if (sym_tag.begin() == sym_tag.end()) {
		sym_tag = root.find_loop("_space_group_symop_operation_xyz");
		sym_string = "_space_group_symop_operation_xyz";
	}

	sym_loop = sym_tag.get_loop();
	sym_pos = sym_loop->find_tag("_symmetry_equiv_pos_as_xyz");

	if (sym_pos == -1)
		sym_pos = 0;

	const int sym_val_size = sym_loop->values.size();
	const int sym_tag_size = sym_loop->tags.size();

	const int sym_size = sym_val_size / sym_tag_size;

	mu::Parser arithmetic_parser;
	
	double dx, dy, dz;

	arithmetic_parser.DefineVar("x", &dx);
	arithmetic_parser.DefineVar("y", &dy);
	arithmetic_parser.DefineVar("z", &dz);

	int num_failures = 0;
	for (int k = 0; k < sym_size; k++) {
		auto sym_group = sym_loop->val(k, sym_pos);
		std::vector<std::string> sym_equations;

		sym_group.erase(std::remove_if(sym_group.begin(), sym_group.end(), boost::is_space()), sym_group.end());

		boost::erase_all(sym_group, "'");
		boost::split(sym_equations, sym_group, boost::is_any_of(","));


		
		for (int i = 0; i < atoms_size; i++) {
			if (!pp.show_hydrogen && atom_loop->val(i, type_pos) == "H")
				continue;

			dx = io::stod(atom_loop->val(i, xyz_pos + 0));
			dy = io::stod(atom_loop->val(i, xyz_pos + 1));
			dz = io::stod(atom_loop->val(i, xyz_pos + 2));

			std::vector<Eigen::Vector3d> temp_atoms;
			const int total_atoms = static_cast<int>(pow(pp.extension, 3));

			Eigen::Vector3d expr_val;
			for (auto j = 0; j < sym_equations.size(); j++) {
				try {
					arithmetic_parser.SetExpr(sym_equations[j]);
					expr_val(j) = arithmetic_parser.Eval();

				}
				catch (mu::Parser::exception_type &e) {
					std::cout << e.GetMsg() << std::endl;
				}
			}


			if (pp.trim_unit_cell) {
				// folding
				for (int j = 0; j < 3; j++) {
					if (expr_val(j) > 1)
						expr_val(j) -= 1;

					if (expr_val(j) < 0)
						expr_val(j) += 1;
				}
			}

			for (int l = 0; l < total_atoms; l++)
				temp_atoms.push_back(expr_val);

			if (pp.extension == 2) {
				// x axis
				temp_atoms[1](0) += 1;
				// y axis
				temp_atoms[2](1) += 1;
				// x-y axis
				temp_atoms[3](0) += 1;
				temp_atoms[3](1) += 1;
				// z axis
				temp_atoms[4](2) += 1;
				// x-z axis
				temp_atoms[5](0) += 1;
				temp_atoms[5](2) += 1;
				// y-z axis
				temp_atoms[6](1) += 1;
				temp_atoms[6](2) += 1;
				// x-y-z axis
				temp_atoms[7](0) += 1;
				temp_atoms[7](1) += 1;
				temp_atoms[7](2) += 1;
			}
			
			std::vector<bool> cf;
			for (int l = 0; l < total_atoms; l++) {
				cf.push_back(false);
			}
			if (pp.trim_unit_cell)
				for (int l = 0; l < total_atoms; l++) {
					for (int j = 0; j < 3; j++)
						if (temp_atoms[l](j) < 0 || temp_atoms[l](j) > pp.extension) {
							cf[l] = true;
						}
				}

			if (pp.trim_unit_cell) {
				for (auto cfv : cf) {
					if (cfv == true)
						num_failures++;
				}
			}

			for (int l = 0; l < total_atoms; l++) {
				if (!cf[l]) {
					cs.atom_names.push_back(atom_loop->val(i, type_pos));
					cs.atom_vectors.push_back(temp_atoms[l]);
				}
			}
		}
	}

	//std::cout << "Num failed trimmings: " << num_failures << std::endl;
}

void perform_orthogonal_conversion(crystal_structure &cs) {
	Eigen::Matrix3d cell = math::fr2cart({ cs.unit_cell.a, cs.unit_cell.b, cs.unit_cell.c, cs.unit_cell.alpha, cs.unit_cell.beta, cs.unit_cell.gamma });

	for (auto &atom : cs.atom_vectors) { atom = cell * atom; }
}


std::vector<Point_3> make_points(crystal_structure& cs, parser_parameters& pp) {
	std::vector<Point_3> lp;

	for (auto atom : cs.atom_vectors) {
		lp.emplace_back(atom(0), atom(1), atom(2));
	}

	return lp;
}

std::vector<double> make_weights(crystal_structure& cs, parser_parameters& pp) {
	std::vector<double> lw;

	for (auto atom : cs.atom_names) {
		lw.push_back(radTable[atom]);
	}

	return lw;
}

void compute_persistence(crystal_structure& cs, parser_parameters& pp, const std::string& off_file_points, const std::string& weight_file, const std::string& output_file_diag, int coeff_field_characteristic, Filtration_value min_persistence, int dimension, double tolerance) {
	// Read the OFF file (input file name given as parameter) and triangulate points
	Gudhi::Points_3D_off_reader<Point_3> off_reader(off_file_points);
	// Check the read operation was correct
	if (!off_reader.is_valid()) {
		std::cerr << "Unable to read OFF file " << off_file_points << std::endl;
		exit(-1);
	}
	// Retrieve the points
	std::vector<Point_3> lp = off_reader.get_point_cloud();

	auto tflp = make_points(cs, pp);//
	auto tflw = make_weights(cs, pp);

	for (int i = 0; i < lp.size(); i++) {
		if (lp[i] != tflp[i]) {
			//lp[i] = tflp[i];
		}
	}

	//lp = tflp;

	// Read weights information from file
	std::ifstream weights_ifstr(weight_file);
	std::vector<Weighted_point_3> wp;
	if (weights_ifstr.good()) {
		double weight = 0.0;
		std::size_t index = 0;
		wp.reserve(lp.size());
		// Attempt read the weight in a double format, return false if it fails
		while ((weights_ifstr >> weight) && (index < lp.size())) {
			wp.push_back(Weighted_point_3(lp[index], tflw[index]));
			index++;
		}
		if (index != lp.size()) {
			std::cerr << "Bad number of weights in file " << weight_file << std::endl;
			exit(-1);
		}
	}
	else {
		std::cerr << "Unable to read weights file " << weight_file << std::endl;
		exit(-1);
	}
	// alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode.
	Alpha_shape_3 as(wp.begin(), wp.end(), 0, Alpha_shape_3::GENERAL);

	// filtration with alpha values from alpha shape
	std::vector<Object> the_objects;
	std::vector<Alpha_value_type> the_alpha_values;
	Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
		std::back_inserter(the_alpha_values));
	as.filtration_with_alpha_values(disp);

	Alpha_shape_3::size_type count_vertices = 0;
	Alpha_shape_3::size_type count_edges = 0;
	Alpha_shape_3::size_type count_facets = 0;
	Alpha_shape_3::size_type count_cells = 0;
	// Loop on objects vector
	Vertex_list vertex_list;
	ST simplex_tree;
	Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
	std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
	std::vector<Weighted_point_3> aspoints;
	for (auto object_iterator : the_objects) {
		// Retrieve Alpha shape vertex list from object
		if (const Cell_handle *cell = CGAL::object_cast<Cell_handle>(&object_iterator)) {
			vertex_list = from_cell<Vertex_list, Cell_handle>(*cell);
			count_cells++;
		}
		else if (const Facet *facet = CGAL::object_cast<Facet>(&object_iterator)) {
			vertex_list = from_facet<Vertex_list, Facet>(*facet);
			count_facets++;
		}
		else if (const Edge_3 *edge = CGAL::object_cast<Edge_3>(&object_iterator)) {
			vertex_list = from_edge<Vertex_list, Edge_3>(*edge);
			count_edges++;
		}
		else if (const Vertex_handle *vertex = CGAL::object_cast<Vertex_handle>(&object_iterator)) {
			count_vertices++;
			vertex_list = from_vertex<Vertex_list, Vertex_handle>(*vertex);
		}

		// Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
		Simplex_tree_vector_vertex the_simplex;
		for (auto the_alpha_shape_vertex : vertex_list) {
			Alpha_shape_simplex_tree_map::iterator the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
			if (the_map_iterator == map_cgal_simplex_tree.end()) {
				// alpha shape not found
				Simplex_tree_vertex vertex = map_cgal_simplex_tree.size();

				the_simplex.push_back(vertex);
				map_cgal_simplex_tree.emplace(the_alpha_shape_vertex, vertex);
			}
			else {
				// alpha shape found
				Simplex_tree_vertex vertex = the_map_iterator->second;
				the_simplex.push_back(vertex);
			}

			aspoints.push_back(the_alpha_shape_vertex->point());
			
		}
		// Construction of the simplex_tree
		Filtration_value filtr = /*std::sqrt*/ (*the_alpha_value_iterator);

		simplex_tree.insert_simplex(the_simplex, filtr);
		if (the_alpha_value_iterator != the_alpha_values.end())
			++the_alpha_value_iterator;
		else
			std::cout << "This shall not happen" << std::endl;
	}

	// Sort the simplices in the order of the filtration
	simplex_tree.initialize_filtration();
	auto stfsr = simplex_tree.filtration_simplex_range();

	//std::cout << "Simplex_tree dim: " << simplex_tree.dimension() << std::endl;
	// Compute the persistence diagram of the complex
	Persistent_cohomology pcoh(simplex_tree, true);
	// initializes the coefficient field for homology
	pcoh.init_coefficients(coeff_field_characteristic);
	pcoh.compute_persistent_cohomology(min_persistence);

	// sort persistence by dimension and value
	cmp_intervals_by_dim_then_length cmp(&simplex_tree);
	auto persistent_pairs = pcoh.get_persistent_pairs();
	std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);
	std::vector<double> hpv;

	for (auto pair : persistent_pairs) {
		//std::cout << simplex_tree.filtration(get<0>(pair)) << " " << simplex_tree.filtration(get<1>(pair)) << std::endl;
		auto svr = simplex_tree.simplex_vertex_range(get<0>(pair));
		std::vector<int> vector_indices;
		for (auto itl : svr) {
			if (itl > -1) vector_indices.push_back(itl);
		}

		double pv = simplex_tree.filtration(get<1>(pair)) - simplex_tree.filtration(get<0>(pair));

		if (simplex_tree.dimension(get<0>(pair)) == dimension)
			hpv.push_back(pv);
	}

	std::vector<double> hpvd;
	double hpd = -1;
	int hpdi = -1;

	hpvd.push_back(0);
	for (int i = 0; i < hpv.size() - 1; i++) {
		hpvd.push_back(hpv[i] - hpv[i + 1]);
		if (hpvd[i] > hpd) {
			hpd = hpvd[i];
			hpdi = i;
		}
	}

	int nchpd = 0;
	std::vector<int> tccp;
	std::vector<int> hpvdix;
	for (int i = 0; i < hpdi; i++) {
		int cccp  = 0;
		while (hpvd[i] < tolerance) {
			i++;
			cccp++;
		}

		if (cccp > 0) {
			tccp.push_back(cccp);
			hpvdix.push_back(i);
		}
			
	}

	std::cout << "+Number of point clusters above widest diagonal gap: " << tccp.size() << ", with pv: " << hpd << std::endl;
	for (int i = 0; i < tccp.size(); i++ ) {
		std::cout << "\t Cluster " << i + 1 <<  " with : " << tccp[i] << " points, pv: " << hpvd[hpvdix[i]] << std::endl;
	}

	// Output the diagram in filediag
	if (output_file_diag.empty()) {
		pcoh.output_diagram();
	}
	else {
		std::cout << "Result in file: " << output_file_diag << std::endl;
		std::ofstream out(output_file_diag + ".xyz");
		std::ofstream persistence(output_file_diag + ".pers");
		//pcoh.output_diagram(out);
		pcoh.output_diagram_xyz(out, persistence, output_file_diag, dimension);
		out.close();
	}
}

void output_off(const std::string& off_file, const std::string& weights_file, crystal_structure& cs, parser_parameters& pp) {
	std::ofstream off(off_file + ".off");
	std::ofstream weights(weights_file + ".weights");

	double minx = 0;
	double miny = 0;
	double minz = 0;

	if (pp.offset_by_min) {
		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			if (cs.atom_vectors.at(i)(0) < minx)
				minx = cs.atom_vectors.at(i)(0);

			if (cs.atom_vectors.at(i)(1) < miny)
				minx = cs.atom_vectors.at(i)(1);

			if (cs.atom_vectors.at(i)(2) < minz)
				minx = cs.atom_vectors.at(i)(2);
		}

		if (minx < miny)
			miny = minx;

		if (minz < miny)
			miny = minz;

		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			cs.atom_vectors.at(i)(0) += -miny;
			cs.atom_vectors.at(i)(1) += -miny;
			cs.atom_vectors.at(i)(2) += -miny;
		}
	}

	if (pp.check_duplicates) {
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;

		std::vector<std::string> tags;

		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			auto atom = cs.atom_vectors.at(i);

			bool flag = true;

			for (auto j = i + 1; j < cs.atom_vectors.size(); j++) {

				auto atom2 = cs.atom_vectors.at(j);

				if (atom(0) == atom2(0) && atom(1) == atom2(1) && atom(2) == atom2(2)) {
					flag = false;
					break;
				}
			}

			if (flag) {
				x.push_back(atom(0));
				y.push_back(atom(1));
				z.push_back(atom(2));

				tags.push_back(cs.atom_names.at(i));
			}
		}

		off << "OFF" << std::endl;
		off << x.size() << " " << 0 << " " << 0 << std::endl;

		for (int i = 0; i < x.size(); i++) {
			off << x[i] << "\t\t" << y[i] << "\t\t" << z[i] << std::endl;
			weights << radTable[tags[i]] << std::endl;
		}

	}
	else {
		off << "OFF" << std::endl;
		off << cs.atom_vectors.size() << " " << 0 << " " << 0 << std::endl;
		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			auto atom = cs.atom_vectors.at(i);
			off << atom(0) << "\t\t" << atom(1) << "\t\t" << atom(2) << std::endl;
			weights << radTable[cs.atom_names[i]] << std::endl;
		}
	}

	
	weights.close();
	off.close();
}

void output_xyz(const std::string& file, crystal_structure& cs, parser_parameters& pp) {
	std::ofstream xyz(file + ".xyz");
	
	double minx = 0;
	double miny = 0;
	double minz = 0;

	if (pp.offset_by_min) {
		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			if (cs.atom_vectors.at(i)(0) < minx)
				minx = cs.atom_vectors.at(i)(0);

			if (cs.atom_vectors.at(i)(1) < miny)
				minx = cs.atom_vectors.at(i)(1);

			if (cs.atom_vectors.at(i)(2) < minz)
				minx = cs.atom_vectors.at(i)(2);
		}

		if (minx < miny)
			miny = minx;

		if (minz < miny)
			miny = minz;

		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			cs.atom_vectors.at(i)(0) += -miny;
			cs.atom_vectors.at(i)(1) += -miny;
			cs.atom_vectors.at(i)(2) += -miny;
		}
	}

	if (pp.check_duplicates) {
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;

		std::vector<std::string> tags;

		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			auto atom = cs.atom_vectors.at(i);

			bool flag = true;

			for (auto j = i+1; j < cs.atom_vectors.size(); j++) {

				auto atom2 = cs.atom_vectors.at(j);

				if (atom(0) == atom2(0) && atom(1) == atom2(1) && atom(2) == atom2(2)) {
					flag = false;
					break;
				}
			}

			if (flag) {
				x.push_back(atom(0));
				y.push_back(atom(1));
				z.push_back(atom(2));

				tags.push_back(cs.atom_names.at(i));
			}
		}

		xyz << x.size() << std::endl;
		xyz << file << std::endl;

		for (int i = 0; i < x.size(); i++) {
			xyz << tags[i] << "\t\t" << x[i] << "\t\t" << y[i] << "\t\t" << z[i] << std::endl;
		}

	} else {
		xyz << cs.atom_vectors.size() << std::endl;
		xyz << file << std::endl;

		for (auto i = 0; i < cs.atom_vectors.size(); i++) {
			auto atom = cs.atom_vectors.at(i);
			xyz << cs.atom_names.at(i) << "\t\t" << atom(0) << "\t\t" << atom(1) << "\t\t" << atom(2) << std::endl;
		}
	}

	xyz.close();
}


/*
void program_options(int argc, char *argv[], std::string &off_file_points, std::string &output_file_diag,
	int &coeff_field_characteristic, Filtration_value &min_persistence) {
	namespace po = boost::program_options;
	po::options_description hidden("Hidden options");
	hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
		"Name of file containing a point set. Format is one point per line:   X1 ... Xd ");
	po::options_description visible("Allowed options", 100);
	visible.add_options()("help,h", "produce help message")(
		"output-file,o", po::value<std::string>(&output_file_diag)->default_value(std::string()),
			"Name of file in which the persistence diagram is written. Default print in std::cout")(
		"field-charac,p", po::value<int>(&coeff_field_characteristic)->default_value(11),
			"Characteristic p of the coefficient field Z/pZ for computing homology.")(
		"min-persistence,m", po::value<Filtration_value>(&min_persistence),
			"Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
			"intervals");
	po::positional_options_description pos;
	pos.add("input-file", 1);
	po::options_description all;
	all.add(visible).add(hidden);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
	po::notify(vm);
	if (vm.count("help") || !vm.count("input-file")) {
		std::cout << std::endl;
		std::cout << "Compute the persistent homology with coefficient field Z/pZ \n";
		std::cout << "of a 3D Alpha complex defined on a set of input points.\n \n";
		std::cout << "The output diagram contains one bar per line, written with the convention: \n";
		std::cout << "   p   dim b d \n";
		std::cout << "where dim is the dimension of the homological feature,\n";
		std::cout << "b and d are respectively the birth and death of the feature and \n";
		std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;
		std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
		std::cout << visible << std::endl;
		exit(-1);
	}
}
*/

int main(int argc, char* argv[])
{
	const std::string cif_path = "data/cif/";
	std::vector<std::string> folder_data;
	
	//folder_data.push_back("ZnAlaPyr_open.cif");
	//folder_data.push_back("ZnAlaPyr_closed.cif");
	folder_data.push_back("ARAHIM.cif");
	folder_data.push_back("DAXNEY.cif");
	folder_data.push_back("DAXNIC.cif");
	folder_data.push_back("FAWCEN.cif"); // -- Crash
	folder_data.push_back("LAJKUE.cif"); // -- Crash
	folder_data.push_back("MESGIB.cif");
	folder_data.push_back("SODDOE.cif");
	folder_data.push_back("WESZOK.cif");
	folder_data.push_back("YOBPOW.cif");
	folder_data.push_back("YOBPUC.cif");
	folder_data.push_back("ZIDYUG.cif");

	initializeRadTable();

	CIFHandler ch;
	crystal_structure cs;
	parser_parameters pp;

	pp.show_hydrogen = true;
	pp.extension = 1;

	pp.trim_unit_cell = true;
	pp.reduce_unit_cell = false;

	pp.check_duplicates = false;
	pp.offset_by_min = false;

	for (std::string file : folder_data) {
		const std::string output_xyz_path = "data/xyz/" + file.substr(0, file.find(".cif"));
		const std::string output_off_path = "data/off/" + file.substr(0, file.find(".cif"));
		const std::string output_weights_path = "data/weights/" + file.substr(0, file.find(".cif"));
		const std::string output_pers_path = "data/pers/" + file.substr(0, file.find(".cif"));

		const std::string xyz_file = output_xyz_path + ".xyz";
		const std::string off_file = output_off_path + ".off";
		const std::string weights_file = output_weights_path + ".weights";
		const std::string persistence_file = output_pers_path + ".pers";

		parse_file(cif_path + file, ch, cs, pp);
		perform_orthogonal_conversion(cs);

		// TODO: use iterator for checking duplicates as it speeds up tremendously
		output_xyz(output_xyz_path, cs, pp);
		output_off(output_off_path, output_weights_path, cs, pp);

		compute_persistence(cs, pp, off_file, weights_file, output_pers_path, 2, 0, 1, 0.01);
	}

	return 0;
}