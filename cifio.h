#ifndef _CIF_IO_H
#define _CIF_IO_H

#include <gemmi/cif.hpp>
#include <gemmi/to_cif.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

namespace cif = gemmi::cif;
using namespace std;

class CIFHandler {
public:


	/*
	*	Copy the CIF Document 
	*	param1: CIF file as Document class
	*	return: gemmi::cif::Document class 
	*/
	cif::Document* copyCIFrom(cif::Document doc);
	
	/*
	*	Change the cell parameters of a CIF Document
	*	param1: Document to modify
	*	param2: vector of 6 parameters in this order { length_a, length_b, length_c, alpha, beta, gamma }
	*	return: gemmi::cif::Document class with parameters
	*/
	bool modifyCIFcellParameters(cif::Document* doc, vector<double> &parameters);

	/*
	*	Perturb the lattice about the quantity specified in parameters_offset
	*	param1: Document to perturb
	*	param2: vector of 6 offsets for present parameters in this order { length_a, length_b, length_c, alpha, beta, gamma }
	*	param3: Characteristic vector of parameters_offset used to mark the elements to use 
	*	return: true if it has been perturbed, false if some pair in the Document has not been found
	*/
	bool perturbTheLattice(cif::Document* doc, vector<double> &parameters_offset, vector<unsigned short> characteristic_vector = { 1,1,1,1,1,1 });

	/*
*	Create a simple CIF Document with cell parameters provided as input
*	param1: vector of 6 parameters in this order { length of a, length of b, length of c, alpha, beta, gamma }
*	return: gemmi::cif::Document Simple class with parameters
*/
	cif::Document* createSimpleCIF(vector<double> &parameters, const char* filename);

private:
	vector<string> cell_parameters_tags{ "_cell_length_a", "_cell_length_b", "_cell_length_c",
							"_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma" };
	string cifexample = "cifexample.cif";

	/*
	*	Create the skeleton of a CIF Document from an existing CIF Document class (only with tags)
	*	param1: CIF file as Document class
	*	return: gemmi::cif::Document class with only tags
	*/
	cif::Document* createEmptyCIFDocumentFrom(cif::Document doc, string filename);

};

#endif // !_CIF_IO_H
