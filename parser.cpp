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
	//unit_cell_strings_.push_back("_cell_volume");

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

	cif::Loop* atom_loop = root.find_loop(atom_site_strings_[0]).get_loop();
	cif::Loop* symm_loop = nullptr;

	std::vector<int> atom_ids;
	for (auto s : atom_site_strings_) {
		atom_ids.push_back(atom_loop->find_tag(s));
	}

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

	const int num_atoms = atom_loop->values.size() / atom_loop->tags.size();
	const int num_symms = symm_loop->values.size() / symm_loop->tags.size();

	mu::Parser mup;

	double dx, dy, dz;

	mup.DefineVar("x", &dx);
	mup.DefineVar("y", &dy);
	mup.DefineVar("z", &dz);

	for (int k = 0; k < num_symms; k++) {
		auto symm_group = symm_loop->val(k, symm_pos);
		std::vector<std::string> symm_eq;

		symm_group.erase(std::remove_if(symm_group.begin(), symm_group.end(), boost::is_space()), symm_group.end());

		boost::erase_all(symm_group, "'");
		boost::split(symm_eq, symm_group, boost::is_any_of(","));

		for (int i = 0; i < num_atoms; i++) {
			if (!p_.show_hydrogen && atom_loop->val(i, atom_ids[0]) == "H")
				continue;

			dx = io::stod(atom_loop->val(i, atom_ids[1]));
			dy = io::stod(atom_loop->val(i, atom_ids[2]));
			dz = io::stod(atom_loop->val(i, atom_ids[3]));

			Eigen::Vector3d expr_val;
			for (int j = 0; j < symm_eq.size(); j++) {
				try {
					mup.SetExpr(symm_eq[j]);
					expr_val(j) = mup.Eval();
				} catch (mu::Parser::exception_type &e) {
					std::cout << e.GetMsg() << std::endl;
				}
			}

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

			const int max_atoms = static_cast<int>(pow(p_.extension, 3));
			std::vector<Eigen::Vector3d> atoms;

			for (int j = 0; j < max_atoms; j++) {
				atoms.push_back(expr_val);
			}

			// manual for now
			if (p_.extension == 2) {
				// x axis
				if (p_.extend_x) {
					atoms[1](0) += 1;
				}
				// y axis
				if (p_.extend_y) {
					atoms[2](1) += 1;
				}
				// x-y axis
				if (p_.extend_x && p_.extend_y) {
					atoms[3](0) += 1;
					atoms[3](1) += 1;
				}
				if (p_.extend_z) {
					// z axis
					if (p_.flatten_z)
						atoms[4](2) += 0.0001;
					else
						atoms[4](2) += 1;
					// z-x axis
					if (p_.extend_x) {
						atoms[5](0) += 1;
						if (p_.flatten_z)
							atoms[5](2) += 0.0001;
						else
							atoms[5](2) += 1;
					}
					// z-y axis
					if (p_.extend_y) {
						atoms[6](1) += 1;
						if (p_.flatten_z)
							atoms[6](2) += 0.0001;
						else
							atoms[6](2) += 1;
					}
					// z-x-y axis
					if (p_.extend_x && p_.extend_y) {
						atoms[7](0) += 1;
						atoms[7](1) += 1;
						if (p_.flatten_z)
							atoms[7](2) += 0.0001;
						else
							atoms[7](2) += 1;
					}	
				}
			}

			std::vector<bool> cf;
			for (int l = 0; l < max_atoms; l++) {
				cf.push_back(false);
			}
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
