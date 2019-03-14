#pragma once

//--- Standard Includes
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

//--- Boost Includes
#include <boost/iostreams/device/mapped_file.hpp>

namespace io {
	static std::string read(const std::string& filepath) {
		boost::iostreams::mapped_file mmap(filepath.c_str(), boost::iostreams::mapped_file::readonly);

		return std::string(mmap.const_data());
	}

	static std::vector<std::string> split(const std::string& file) {
		std::istringstream input(file);
		std::vector<std::string> temp;

		for(std::string l; std::getline(input, l, '\n');) {
			temp.push_back(l);
		}

		return temp;
	}

	static std::string read_sequential(const std::string& filepath) {
		std::ifstream infile{ filepath };

		infile.open(filepath);

		std::string file_contents{ std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>() };

		infile.close();

		return file_contents;
	}

	static std::vector<std::string> tokenize(const std::string line, const std::string& delim) {
		std::vector<std::string> token;
		auto temp = line;

		while (!temp.empty()) {
			int ndx = static_cast<int>(temp.find_first_of(delim));
			if (ndx > 0) {
				token.push_back(temp.substr(0, ndx));
			}
			else if (ndx == -1) {
				token.push_back(temp);
				return token;
			}

			temp = temp.substr(ndx + 1);
		}

		return token;
	}

	static double stod(std::string str) {
		std::istringstream i(str);
		double x;
		
		if (!(i >> x)) {
			std::cout << "Bad string to double conversion" << std::endl;
			exit(0);
		}

		return x;
	}
}
