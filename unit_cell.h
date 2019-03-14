#pragma once

#include <vector>
#include <string>

#include "math.h"
#include "file.h"

class unit_cell
{
public:
	double a, b, c;
	double alpha, beta, gamma;
	double volume;

	double A, B, C, xi, eta, zeta;
	unit_cell() {
		a = b = c = alpha = beta = gamma = volume = 0;
	}

	unit_cell(double a, double b, double c, double alpha, double beta, double gamma) {
		this->a = a;
		this->b = b;
		this->c = c;
		this->alpha = alpha;
		this->beta = beta;
		this->gamma = gamma;
		volume = 0;
	}

	void init(std::vector<std::string> file) {
		if(!file.empty()) {
			load_file(file);
		}
		else {
			a = b = c = alpha = beta = gamma = volume = 0;
		}
	}

	void prepare_reduction() {
		A = pow(a, 2);
		B = pow(b, 2);
		C = pow(c, 2);

		//boost::math::cos

		xi = 2 * b * c * boost::math::cos_pi(alpha / 180.f);
		eta = 2 * a * c * boost::math::cos_pi(beta / 180.f);
		zeta = 2 * a * b * boost::math::cos_pi(gamma / 180.f);
	}

private:
	void load_file(std::vector<std::string> file) {
		std::vector<std::string> list;
		list.emplace_back("data_");
		list.emplace_back("_cell_length_a");
		list.emplace_back("_cell_length_b");
		list.emplace_back("_cell_length_c");
		list.emplace_back("_cell_angle_alpha");
		list.emplace_back("_cell_angle_beta");
		list.emplace_back("_cell_angle_gamma");
		list.emplace_back("_cell_volume");

		for (const auto& l : file) {
			auto token = io::tokenize(l, " ()\r\t");

			if (!token.empty()) {
				if (token[0] == list[1])
					a = io::stod(token[1]);
				else if (token[0] == list[2])
					b = io::stod(token[1]);
				else if (token[0] == list[3])
					c = io::stod(token[1]);
				else if (token[0] == list[4])
					alpha = io::stod(token[1]);
				else if (token[0] == list[5])
					beta = io::stod(token[1]);
				else if (token[0] == list[6])
					gamma = io::stod(token[1]);
				else if (token[0] == list[7])
					volume = io::stod(token[1]);
			}
		}
	}
};
