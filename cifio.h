#pragma once

#include <gemmi/cif.hpp>
#include <gemmi/to_cif.hpp>

namespace cif = gemmi::cif;

class CIFHandler {
public:

	/*
	*	Copy the CIF Document 
	*	param1: CIF file as Document class
	*	return: gemmi::cif::Document class 
	*/
	cif::Document* copyCIFrom(cif::Document doc);

private:
	std::vector<std::string> cell_parameters_tags{ "_cell_length_a", "_cell_length_b", "_cell_length_c",
							"_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma" };
};