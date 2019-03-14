#include "cifio.h"

cif::Document* CIFHandler::createEmptyCIFDocumentFrom(cif::Document doc, string filename) {

	cif::Document* new_doc = new cif::Document();
	new_doc->source = filename;
	vector<string> curr_row;
	int spec_block_number = 0;
	for (cif::Block &block : doc.blocks) {

		cif::Block* new_block = (block.name.compare("global") == 0) ? &(new_doc->add_new_block(block.name)) : &(new_doc->add_new_block(string("specifics_") + to_string(spec_block_number++)));

		for (cif::Item& item : block.items) {
			if (item.type == cif::ItemType::Pair) {
				new_block->set_pair(item.pair.at(0), string("\t?"));
			}
			if (item.type == cif::ItemType::Loop) {
				cif::Loop* new_loop;
				int tags_number = item.loop.tags.size();
				new_loop = &(new_block->init_loop(string(), item.loop.tags));
				// At least one add_row must be used to print the tags
				for (int j = 0; j < tags_number; ++j) {
					curr_row.push_back(string("?"));
				}
				new_loop->add_row(curr_row);
				curr_row.clear();
			}
		}
	}
	return new_doc;
}

cif::Document* CIFHandler::copyCIFrom(cif::Document doc) {

	cif::Document* new_doc = new cif::Document();
	new_doc->source = doc.source;
	vector<string> curr_row;
	for (cif::Block &block : doc.blocks) {
		cif::Block* new_block = &(new_doc->add_new_block(block.name));

		for (cif::Item& item : block.items) {
			if (item.type == cif::ItemType::Pair) {
				new_block->set_pair(item.pair.at(0), item.pair.at(1));
			}
			if (item.type == cif::ItemType::Loop) {
				cif::Loop* new_loop;
				int tags_number = item.loop.tags.size(), row_number = item.loop.values.size() / item.loop.tags.size();
				new_loop = &(new_block->init_loop(string(), item.loop.tags));
				
				// This part copies values from a doc to the new one
				for (int i = 0; i < row_number; ++i) {
					for (int j = 0; j < tags_number; ++j) {
						curr_row.push_back(item.loop.values[i*tags_number+j]);
					}
					new_loop->add_row(curr_row);
					curr_row.clear();
				}				
			}
		}
	}
	return new_doc;
}


bool CIFHandler::modifyCIFcellParameters(cif::Document* doc, vector<double> &parameters) {

	if (parameters.size() != 6) {
		cerr << "Error in modifyCIF: The vector must have 6 parameters, 3 lengths and 3 angles" << endl;
		return false;
	}
	bool found = false;
	for (cif::Block &block : doc->blocks) {
		for (unsigned short i = 0; i < this->cell_parameters_tags.size(); ++i) {
			const cif::Pair* p = block.find_pair(this->cell_parameters_tags[i]);
			if (p != nullptr) {
				block.set_pair(this->cell_parameters_tags[i], to_string(parameters[i]));
				found = true;
			} else {
				// no pair with tags
				found = false;
				break;
			}
		}
	}
	if (!found) {
		cerr << "No pairs have been found" << endl;
	}
	return found;
}


bool CIFHandler::perturbTheLattice(cif::Document* doc, vector<double> &parameters_offset, vector<unsigned short> characteristic_vector) {
	if ( (parameters_offset.size() != 6 ) && ( characteristic_vector.size() != 6 ) ) {
		cerr << "Error in perturbTheLattice: The vectors must have 6 parameters" << endl;
		return false;
	}
	bool found = false;
	for (cif::Block &block : doc->blocks) {
		for (unsigned short i = 0; i < this->cell_parameters_tags.size(); ++i) {
			const cif::Pair* p = block.find_pair(this->cell_parameters_tags[i]);
			if (p != nullptr) {
				if (characteristic_vector[i] == 1) {
					block.set_pair(this->cell_parameters_tags[i], to_string( stod( p->at(1) ) + parameters_offset[i]));
				}
				found = true;
			}
			else {
				// no pair with tags
				found = false;
				break;
			}
		}
	}
	if (!found) {
		cerr << "No pairs have been found" << endl;
	}
	return found;
}

cif::Document* CIFHandler::createSimpleCIF(vector<double> &parameters, const char* filename) {

	cif::Document* new_doc, doc_example;

	if (parameters.size() != 6) {
		cerr << "Error in createSimpleCIF: The vector must have 6 parameters, 3 lengths and 3 angles" << endl;
		return nullptr;
	}

	doc_example = cif::read_file(this->cifexample);
	new_doc = this->createEmptyCIFDocumentFrom(doc_example, string(filename));
	this->modifyCIFcellParameters(new_doc, parameters);

	return new_doc;

}
