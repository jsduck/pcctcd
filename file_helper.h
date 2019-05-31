#pragma once

//--- Standard Includes
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

//--- Boost Includes
#include <boost/iostreams/device/mapped_file.hpp>

namespace io {
	/*
	 * Read a file and load it into heap directly
	 * param1: string of the file path
	 * return: string of the full content of the file
	 */
	static std::string read(const std::string& filepath) {
		boost::iostreams::mapped_file mmap(filepath.c_str(), boost::iostreams::mapped_file::readonly);

		return std::string(mmap.const_data());
	}

	/*
	 * Split a string based on \newlines
	 * param1: string with desired data to be handled
	 * return: vector containing the separated contents of the input
	 */
	static std::vector<std::string> split(const std::string& file) {
		std::istringstream input(file);
		std::vector<std::string> temp;

		for(std::string l; std::getline(input, l, '\n');) {
			temp.push_back(l);
		}

		return temp;
	}

	/*
	 * Read a an entire file sequentially
	 * param1: string of the file path
	 * return: string of the entire file contents
	 */
	static std::string read_sequential(const std::string& filepath) {
		std::ifstream infile{ filepath };

		infile.open(filepath);

		std::string file_contents{ std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>() };

		infile.close();

		return file_contents;
	}

	/*
	 * Split a string based on given delimitators
	 * param1: string to be tokenised
	 * param2: string containging delimitators e.g. "?!@"
	 * return: std vector containing the tokenized input
	 */
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

	/*
	 * Conversion of string to double value using streams
	 * param1: string to be converted
	 * return: double value
	 */
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
