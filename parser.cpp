//-- Base Include
#include "parser.h"

//-- Local Includes
#include "cifio.h"

//-- Helper Includes
#include "file_helper.h"

//-- muParser Includes
#include <muParser.h>
#include <algorithm/string/classification.hpp>

//-- Boost Includes
#include <boost/version.hpp>
#include <boost/program_options.hpp>
#include <boost/variant.hpp>
#include <boost/algorithm/string.hpp>

#include <Eigen/Eigen>

Parser::Parser() : p_{} {
	unit_cell_strings_.emplace_back("_cell_length_a");
	unit_cell_strings_.emplace_back("_cell_length_b");
	unit_cell_strings_.emplace_back("_cell_length_c");
	unit_cell_strings_.emplace_back("_cell_angle_alpha");
	unit_cell_strings_.emplace_back("_cell_angle_beta");
	unit_cell_strings_.emplace_back("_cell_angle_gamma");
	unit_cell_strings_.emplace_back("_cell_volume");

	//atom_site_strings_.emplace_back("_atom_site_label");
	atom_site_strings_.emplace_back("_atom_site_type_symbol");
	atom_site_strings_.emplace_back("_atom_site_fract_x");
	atom_site_strings_.emplace_back("_atom_site_fract_y");
	atom_site_strings_.emplace_back("_atom_site_fract_z");

	symm_loop_strings_.emplace_back("_symmetry_equiv_pos_as_xyz");
	symm_loop_strings_.emplace_back("_space_group_symop_operation_xyz");
}

void Parser::set_options(options& p) {
	p_ = p;
}

//template<typename T>
void Parser::read(std::string filepath, Crystal& t) const {
	CIFHandler cif;
	// Cell extraction and ops
	auto root = cif.copyCIFrom(cif::read_file(filepath))->blocks[0];
	{
		std::vector<double> v;
		for (auto s : unit_cell_strings_) {
			v.push_back(io::stod(*root.find_value(s)));
		}

		t.append_cell(v);
	}
	if (p_.reduce_unit_cell) {
		t.reduce_cell(0.0001);
	}

	// Extract loop data
	cif::Loop* atom_loop = root.find_loop(atom_site_strings_[0]).get_loop();
	cif::Loop* symm_loop = nullptr;

	// Exctrat atom data
	std::vector<int> atom_ids;
	for (auto s : atom_site_strings_) {
		atom_ids.push_back(atom_loop->find_tag(s));
	}

	// As 2 symm loops have been found in test data
	// this part tests which one is present in the CIF file
	int symm_pos = 0;
	for (auto s : symm_loop_strings_) {
		cif::Column symm_col = root.find_loop(s);

		if (symm_col.begin() != symm_col.end()) {
			symm_loop = symm_col.get_loop();

			symm_pos = symm_loop->find_tag(s);
			if (symm_pos == -1)
				symm_pos = 0;

			break;
		}
	}

	// Expected amount of atoms and symmetry ops
	const int num_atoms = atom_loop->values.size() / atom_loop->tags.size();
	const int num_symms = symm_loop->values.size() / symm_loop->tags.size();

	mu::Parser mup;
	double dx, dy, dz;

	mup.DefineVar("x", &dx);
	mup.DefineVar("y", &dy);
	mup.DefineVar("z", &dz);

	for (int k = 0; k < num_symms; k++) {
		// Exact symmetry op data and erase all undesired noise such as whitespace or punctuation
		auto symm_group = symm_loop->val(k, symm_pos);
		std::vector<std::string> symm_eq;

		symm_group.erase(std::remove_if(symm_group.begin(), symm_group.end(), boost::is_space()), symm_group.end());

		boost::erase_all(symm_group, "'");
		boost::split(symm_eq, symm_group, boost::is_any_of(","));

		for (int i = 0; i < num_atoms; i++) {
			// Skip hydrogens if switch is enabled at runtime
			if (!p_.show_hydrogen && atom_loop->val(i, atom_ids[0]) == "H")
				continue;

			// Extract xyz fractional coordss
			dx = io::stod(atom_loop->val(i, atom_ids[1]));
			dy = io::stod(atom_loop->val(i, atom_ids[2]));
			dz = io::stod(atom_loop->val(i, atom_ids[3]));

			// Interpret symmetry strings
			Eigen::Vector3d expr_val;
			for (int j = 0; j < symm_eq.size(); j++) {
				try {
					mup.SetExpr(symm_eq[j]);
					expr_val(j) = mup.Eval();
				} catch (mu::Parser::exception_type &e) {
					std::cout << e.GetMsg() << std::endl;
				}
			}

			// Move outside atoms outside of 0-1 range back in the unit cell
			if (p_.trim_unit_cell) {
				for (int j = 0; j < 3; j++) {
					if (expr_val(j) > 1)
						expr_val(j) -= 1;

					if (expr_val(j) < 0)
						expr_val(j) += 1;
				}

				if (p_.flatten_z)
					expr_val(2) = 0;
			}

			// Create temporary storage for atoms based on desired dimension
			const int max_atoms = static_cast<int>(pow(p_.extension, 3));
			std::vector<Eigen::Vector3d> atoms;

			for (int j = 0; j < max_atoms; j++) {
				atoms.push_back(expr_val);
			}

			// Extension based on dimension used as input, support any number > 0
			for (int z_i = 0; z_i < p_.extension; z_i++) {
				for (int y_i = 0; y_i < p_.extension; y_i++) {
					for (int x_i = 0; x_i < p_.extension; x_i++) {
						auto index = pow(p_.extension, 2) * z_i + p_.extension * y_i + x_i;

						atoms[index](0) += x_i;
						atoms[index](1) += y_i;
						atoms[index](2) += z_i;
					}
				}
			}

			std::vector<bool> cf;
			for (int l = 0; l < max_atoms; l++) {
				cf.push_back(false);
			}
			// If atom is position outside estimated boundries do not add it in the point cloud
			if (p_.trim_unit_cell)
				for (int l = 0; l < max_atoms; l++) {
					for (int j = 0; j < 3; j++)
						if (atoms[l](j) < 0 || atoms[l](j) > p_.extension) {
							cf[l] = true;
						}
				}

			std::string tag = atom_loop->val(i, atom_ids[0]);
			for (int j = 0; j < max_atoms; j++) {
				if (!cf[j]) t.append_atom(tag, atoms[j]);
			}
		}
	}
}
