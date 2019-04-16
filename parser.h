#pragma once

//-- Standard Includes
#include <string>
#include <vector>

//-- Boost Includes

//-- Local Includes
#include "crystal_str.h"

class Parser
{
public:
	struct options
	{
		bool show_hydrogen,
			trim_unit_cell,
			reduce_unit_cell,
			check_duplicates,
			extend_x,
			extend_y,
			extend_z,
			flatten_z;
		unsigned int extension;
	};
public:
	Parser();

	//template<typename T>
	void read(std::string filepath, Crystal& t) const;

	void set_options(options& p);

private:
	options p_;

	std::vector<std::string> unit_cell_strings_;
	std::vector<std::string> atom_site_strings_;
	std::vector<std::string> symm_loop_strings_;
};