//-- Standard Includes
#include <vector>
#include <filesystem>

//-- Local Includes
#include "atom_data.h"
#include "parser.h"
#include "crystal_str.h"
#include "persistence.h"

int main(int argc, char* argv[])
{
	const std::string cif_path(argc >= 2 ? argv[1] : "data/");
	std::vector<std::string> folder_data;

	for (const auto& f : std::experimental::filesystem::directory_iterator(cif_path)) {
		std::string s = f.path().string();

		folder_data.push_back(s.substr(s.find('\\') + 1));
	}

	atom::initializeRadTable();

	Parser::options op{};

	op.show_hydrogen = false;
	op.extension = 2;

	op.trim_unit_cell = true;
	op.reduce_unit_cell = false;

	op.check_duplicates = false;

	op.extend_x = true;
	op.extend_y = true;
	op.extend_z = true;

	op.flatten_z = false;

	Parser p;
	p.set_options(op);

	//for(int i = 0; i < 10;i++)
	std::ofstream reset_test("result_main.out");
	reset_test.close();

	for (std::string file : folder_data) {
		//const std::string xyz_file = "out\\" + file + ".xyz";
		//const std::string off_file = "out\\" + file + ".off";
		//const std::string weights_file = "out\\" + file + ".weights";
		//const std::string persistence_file = "out\\" + file + ".pers";

		Crystal c;
		p.read(cif_path + file, c);

		c.orthogonal();

		//c.output_xyz(xyz_file);
		//c.output_off(off_file);

		Persistence per(c);

		per.calculate("", "", 2, 0);
		//per.calculate_clusters(0.8);

		per.display_complex();

		//std::cout << "-Finished computing -- " << file << std::endl;
		//std::ofstream f;
		//f.open("result_main.out", std::ofstream::out | std::ofstream::app);
		//f << file << "\t\t";
		//f.close();
		//per.output("result_main.out");
		//auto pairs = per.get_pairs();
	}

	return 0;
}