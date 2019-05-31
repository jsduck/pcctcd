#include "cifio.h"

cif::Document* CIFHandler::copyCIFrom(cif::Document doc) {

	cif::Document* new_doc = new cif::Document();
	new_doc->source = doc.source;
	std::vector<std::string> curr_row;
	for (cif::Block &block : doc.blocks) {
		cif::Block* new_block = &(new_doc->add_new_block(block.name));

		for (cif::Item& item : block.items) {
			if (item.type == cif::ItemType::Pair) {
				new_block->set_pair(item.pair.at(0), item.pair.at(1));
			}
			if (item.type == cif::ItemType::Loop) {
				cif::Loop* new_loop;
				int tags_number = item.loop.tags.size(), row_number = item.loop.values.size() / item.loop.tags.size();
				new_loop = &(new_block->init_loop(std::string(), item.loop.tags));
				
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