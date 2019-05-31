//-- Standard Includes
#include <vector>
#include <filesystem>

//-- Local Includes
#include "atom_data.h"
#include "parser.h"
#include "crystal_str.h"
#include "persistence.h"

#include "gudhi/Persistence_intervals.h"
#include "niggli.h"

//-- Boost Includes
#include <boost/filesystem.hpp>
//#define COMPUTE_ALL
//#define ENABLE_TIMER
//#ifdef ENABLE_TIMER
#include <boost/timer/timer.hpp>
using namespace boost::timer;
//#endif
#include <boost/program_options.hpp>

int main(int argc, char* argv[])
{
	atom::initializeRadTable();

	Parser::options op{};
	using namespace boost::program_options;
	options_description desc("Commands");
	desc.add_options()
		("help", "produce command list")
		("input", value<std::string>(),"input path to data")
		("output", value<std::string>(), "output path to data")
		("d", value<int>(), "space used for computing the extended unit cell")
		("h", "option for the exclusion of hydrogen atoms")
		("n", "option for using the non-reduced cell")
		("cl", "option for calculation clusters")
		("dia", "output gnuplot diagram")
		("vis", "interactive visualisation of the alpha complex")
		("tetra", "discard input and use a simple tetrahedron for testing");

	variables_map vm;
	store(parse_command_line(argc, argv, desc), vm);
	notify(vm);

	//std::vector<double> sk{2.00, 20.025, 21.000, 0.05, 90.00, 90.00};
	//std::cout << sk;
	//niggli::reduce(sk, 0.0001);
	//std::cout << sk[0] << " | " << sk[1] << " | " << sk[2] << " | " << sk[3] << " | " << sk[4] << " | " << sk[5] << " | ";

	//-- Setup initial options values
	op.show_hydrogen = true;
	op.extension = 1;

	op.trim_unit_cell = true;
	op.reduce_unit_cell = true;

	op.check_duplicates = false;

	op.extend_x = true;
	op.extend_y = true;
	op.extend_z = true;

	op.flatten_z = false;

	bool single = false;
	bool dia = false;
	bool vis = false;
	bool cls = false;
	bool tetra = false;
	std::string input;
	std::string output;
	if (vm.count("help")) {
		std::cout << desc;
		return 1;
	}
	if (vm.count("input")) {
		input = vm["input"].as<std::string>();
		//std::cout << input.find(".cif") << std::endl;
		single = (input.find(".cif") != std::string::npos);
	} else {
		single = true;
		input = "";
	}
	if (vm.count("output")) {
		output = vm["output"].as<std::string>();
	} else {
		output = "out/";
	}
	if (vm.count("d")) {
		//std::cout << vm["d"].as<int>() << std::endl;
		op.extension = vm["d"].as<int>();
	}
	if (vm.count("h")) {
		op.show_hydrogen = false;
	}
	if (vm.count("n")) {
		op.reduce_unit_cell = false;
	}
	if (vm.count("dia")) {
		dia = true;
	}
	if (vm.count("cl")) {
		cls = true;
	}
	if (vm.count("vis") && single) {
		vis = true;
	}
	if (vm.count("tetra")) {
		tetra = true;
	}

	std::vector <std::string> folder_data;
	std::vector <std::string> folder_path;
	//std::cout << "Single toggle: " << single << std::endl;
	//-- Testing block for manual file input, disabled by default
	bool enable_debug = false;
	if (enable_debug) {
		tetra = false;
		dia = false;
		cls = false;
		vis = true;
		op.show_hydrogen = false;
		op.reduce_unit_cell = false;
		op.extension = 2;

		single = true;
		input = "data/ZnAlaPyr_open.cif";
	}

	//-- Change file loading behaviour based on the extension of the input
	if (single) {
		folder_data.push_back(input.substr(input.find_last_of('\\') + 1));
		folder_path.push_back(input);
	} else {
		for (const auto& f : std::experimental::filesystem::directory_iterator(input)) {
			std::string s = f.path().string();
			folder_data.push_back(s.substr(s.find_last_of('\\') + 1));
			folder_path.push_back(s);
		}

	}

	
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Number of files: " << folder_data.size() << std::endl;
	std::cout << "Extended space: " << op.extension << std::endl;
	std::cout << "Reduced unit cell: " << (op.reduce_unit_cell) << std::endl;
	std::cout << "Calculation clusters: " << (cls) << std::endl;
	std::cout << "Show hydrogens: " << (op.show_hydrogen) << std::endl;
	std::cout << "Show diagram: " << (dia) << std::endl;
	std::cout << "Show visualisation: " << (vis) << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;

	Parser p;
	p.set_options(op);

	//for(int i = 0; i < 10;i++)
	std::ofstream reset_test("results_overview.out");
	reset_test.close();

	bool enable_commands = false;
	int loop_factor = 1;
	int loop_iterator = 0;

	//-- Only used for manually exporting data during testing
	//-- disabled by default
	bool data_export = false;
	std::vector<std::string> data_vector;

	double factor = 3.f;

	cpu_timer timer;
	//-- The do-while loop exists only for testing and is used when
	//-- wanting to run multiple test runs during the same session
	//-- The loop_factor variable declared above controls the amount
	//-- of runs
	do {
		for (int i = 0; i < folder_data.size(); i++) {
			//cpu_timer timer;
			std::string file = folder_data[i];
			std::string path = folder_path[i];
			//const std::string xyz_file = "out\\" + file + ".xyz";
			const std::string off_file = "out\\" + file + ".off";
			//const std::string weights_file = "out\\" + file + ".weights";
			//const std::string persistence_file = "out\\" + file + ".pers";

			Crystal c;
			//std::cout << "-Atempting to read " << path << std::endl;
			//-- Read the file and convert it as planned if
			//-- tetra flag is not set, otherwise create a
			//-- dummy tetrahedron of a given factor
			if (!tetra) {
				p.read(path, c);
				c.orthogonal();
			}
			else {
				c.append_cell({ 1, 1, 1, 90, 90, 90 });
				c.append_atom("null", Eigen::Vector3d(factor * 0, factor * 0, factor * 0), 0);
				c.append_atom("null", Eigen::Vector3d(factor * 2, factor * 0, factor * 0), 0);
				c.append_atom("null", Eigen::Vector3d(factor * 1, factor*sqrt(3), factor * 0), 0);
				c.append_atom("null", Eigen::Vector3d(factor * 1, factor*sqrt(3) / 3.f, factor*sqrt(8.f / 3)), 0);
				c.make_points();
				//std::cout << c.get_atoms().size() << std::endl;
			}


			//c.output_xyz(xyz_file);
			//c.output_off(off_file);

			Persistence per(c);

			//-- First 2 params are emtpy as they are deprecated
			per.calculate("", "", 2, 0);

			//-- Currently NOT IN USE, needs fixing
			if (cls) {
				per.calculate_clusters(0.8);
			}

			if (dia) {
				Gudhi::Persistence_representations::Persistence_intervals intervals(per.get_pairs());
				intervals.plot(file.c_str());
			}
			if (vis) {
				per.display_complex();
			}

			std::cout << "-Finished computing -- " << file << std::endl;

			//-- Create directories if not already existing
			//-- and output pairs data
			std::ofstream f, o;
			boost::filesystem::path dir(output.c_str());
			boost::filesystem::create_directory(dir);
			per.output_pairs(output + file.substr(0, file.find(".cif")) + ".pairs");
			
			f.open("results_overview.out", std::ofstream::out | std::ofstream::app);
			f << file << "\t\t";
			f.close();
			per.output("results_overview.out");
			//-- While the clean method good for memory management, its terrible for
			//-- speed, so it is disabled until a better method is found
			//per.clean();
		}

		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "-Execution completed in: " << timer.format()/* / 1000000000.f << " seconds" */ << std::endl;

		if (data_export) {
			data_vector.push_back(timer.format());
		}

		loop_iterator++;
	} while (loop_iterator < loop_factor);
	boost::filesystem::remove("temp");

	if (data_export) {
		std::ofstream de("data_export");
		for (auto d : data_vector) {
			de << d << std::endl;
		}
		de.close();
	}

	return 0;
}