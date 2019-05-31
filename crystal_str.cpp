//-- Base Include
#include "crystal_str.h"

//-- Local Includes
#include "niggli.h"
#include "atom_data.h"
#include "math_helper.h"

#include <fstream>

//#define TEST0 
#ifdef TEST0
#define FACTOR 3.f
#endif

Crystal::Crystal()
{
	clear();
}

void Crystal::append_cell(std::vector<double> v) {
	for (auto& d : v) {
		unit_cell_.push_back(d);
	}
}

void Crystal::reduce_cell(double epsilon, int loops) {
	niggli::reduce(unit_cell_, epsilon, loops);
}

std::vector<double>& Crystal::get_cell() {
	return unit_cell_;
}

void Crystal::append_atom(std::string label, Eigen::Vector3d coords, double weight)
{
	atoms_.emplace_back(label, coords);
	if (weight != -1)
		weights_.push_back(weight);
	else
		weights_.push_back(pow(atom::radTable[label], 2));
}

std::pair<std::string, Eigen::Vector3d>& Crystal::get_atom(int index)
{
	return atoms_[index];
}

std::vector<std::pair<std::string, Eigen::Vector3d>>& Crystal::get_atoms()
{
	return atoms_;
}

void Crystal::orthogonal() {
	const Eigen::Matrix3d cell = math::fr2cart(unit_cell_);
	points_.clear();

	for (auto& atom : atoms_) {
		atom.second = cell * atom.second;

		points_.emplace(Point_3(atom.second(0), atom.second(1), atom.second(2)), pow(atom::radTable[atom.first], 2));
	}
}

std::vector<Point_3> Crystal::get_points() {
	std::vector<Point_3> t;

#ifdef TEST0
	//t.emplace_back(-1, -1, -1);
	t.emplace_back(FACTOR*0, FACTOR * 0, FACTOR * 0);
	t.emplace_back(FACTOR * 2, FACTOR * 0, FACTOR * 0);
	t.emplace_back(FACTOR * 1, FACTOR*sqrt(3), FACTOR * 0);
	t.emplace_back(FACTOR * 1, FACTOR*sqrt(3)/3.f, FACTOR*sqrt(8.f/3));

#else
	for (auto& atom : atoms_) {
		t.emplace_back(atom.second(0), atom.second(1), atom.second(2));
	}
#endif
	
	return t;
}

std::map<Point_3, double> Crystal::get_orthogonal_points() const {
	return points_;
}

std::vector<double> Crystal::get_vdW_weights() const {
	//std::vector<double> w;
#ifdef TEST0
	//w.push_back(4.f);
	//w.push_back(4.0f);
	//w.push_back(4.0f);
	//w.push_back(0.0f);
#else
	//double max = 0;
	//for (auto& atom : atoms_) {
	//	if (atom::radTable[atom.first] > max) {
	//		max = atom::radTable[atom.first];
	//	}
	//	w.push_back(pow(atom::radTable[atom.first], 2));
	//}
	//std::cout << max;
#endif

	return weights_;
}

void Crystal::output_off(std::string filename)
{
	std::ofstream off(filename);
	off << "OFF" << std::endl;
#ifdef TEST0
	off << 4 << " " << 0 << " " << 0 << std::endl; 
	//off << -1 << " " << -1 << " " << -1 << std::endl;
	off << FACTOR * 0 << " " << FACTOR * 0 << " " << FACTOR * 0 << std::endl;
	off << FACTOR * 2 << " " << FACTOR * 0 << " " << FACTOR * 0 << std::endl;
	off << FACTOR * 1 << " " << FACTOR * sqrt(3) << " " << FACTOR * 0 << std::endl;
	off << FACTOR * 1 << " " << FACTOR * sqrt(3)/3.f << " " << FACTOR * sqrt(8.f/3) << std::endl;
#else
	off << atoms_.size() << " " << 0 << " " << 0 << std::endl;
	//off << points_.size() << " " << 0 << " " << 0 << std::endl;

	for (auto atom : atoms_) {
		//off << std::setprecision(6) << std::fixed << atom.second(0) << "\t\t" << atom.second(1) << "\t\t" << atom.second(2) << std::endl;
		off << atom.second(0) << "\t\t" << atom.second(1) << "\t\t" << atom.second(2) << std::endl;
	}
	for (auto p : points_) {
		//off << p.first[0] << "\t\t" << p.first[1] << "\t\t" << p.first[2] << std::endl;
	}
#endif

	off.close();
}

void Crystal::output_xyz(std::string filename)
{
	std::ofstream xyz(filename);

	xyz << atoms_.size() << std::endl;
	xyz << filename << std::endl;

	for (auto atom : atoms_) {
		xyz << atom.first << "\t\t" << std::setprecision(6) << std::fixed << atom.second(0) << "\t\t" << atom.second(1) << "\t\t" << atom.second(2) << std::endl;
	}

	xyz.close();
}
