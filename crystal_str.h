#pragma once

//-- Standard Includes
#include <vector>
#include <set>

//-- Eigen Includes
#include <Eigen/Eigen>

//-- Local Includes
#include "cgal_helper.h"

class Crystal
{
public:
	/*
	 * Default constructor
	 */
	Crystal();

	/*
	 * Empty all member arrays
	 */
	void clear() {
		unit_cell_.clear();
		atoms_.clear();
		points_.clear();
		weights_.clear();
	}

	/*
	 * Update the unit cell vector
	 * param1: array of doubles representing unit cell values
	 */
	void append_cell(std::vector<double> v);
	/*
	 * Reduce the unit cell using niggli's reduction
	 * param1: tolerance -- DEPRECATED --
	 * param2: max attempts at reducing the unit cell, defaults to 50
	 */
	void reduce_cell(double epsilon, int loops = 50);

	/*
	 * Retrieve the unit cell values
	 * return: array of doubles represeting the unit cell
	 */
	std::vector<double>& get_cell();

	/*
	 * Add an atom to the atom array
	 * param1: a string object represting its atomic symbol
	 * param2: Vector3d of its coordinates
	 * param3: double value representing weight--DEPRECATED-- only used for testing
	 */
	void append_atom(std::string label, Eigen::Vector3d coords, double weight = -1);

	/*
	 * Retrieve an atom from the array
	 * param1: index of the desired atom
	 * return: std::pair of the atom's Vector of position and atomic symbol
	 */
	std::pair<std::string, Eigen::Vector3d>& get_atom(int index);
	/*
	 * Retrieve the vector of atoms
	 * return: vector containing all atoms
	 */
	std::vector<std::pair<std::string, Eigen::Vector3d>>& get_atoms();

	/*
	 * Transform all atoms from fractional to orthogonal coordinates
	 */
	void orthogonal();

	/*
	 * Forces initialisation of the point_map used in visualisation
	 * Under normal circumstances the point_map is created when converting
	 * coordinates from fractional to orthogonal and as such, this method
	 * should not be used unless using user defined custom atom positioning
	 * not loaded form a fractional unit cell
	 */
	inline void make_points() {
		points_.clear();
		point_set_.clear();

		//std::cout << "Point size preortho: " << points_.size() << std::endl;
		int i = 0;
		for (auto& atom : atoms_) {
			points_.emplace(Point_3(atom.second(0), atom.second(1), atom.second(2)), weights_[i++]);
		}
	}
	/*
	 * Retrieve a vector of CGAL Point_3 objects from the point_map
	 * return: vector of Point_3
	 */
	std::vector<Point_3> get_points();
	/*
	 * Retrieve the point_map of orthogonal points
	 * return: map of Point_3 and double, representing weighted points
	 */
	std::map<Point_3, double> get_orthogonal_points() const;

	/*
	 * Retrieve the van der waals weights corresponding to the atom vector
	 * return: vector of weights
	 */
	std::vector<double> get_vdW_weights() const;

	/*
	 * Output point map for a .off extension format
	 */
	void output_off(std::string filename);
	/*
	 * Output point map for a .xyz extension format
	 */
	void output_xyz(std::string filename);

	/* 
	 * Retrieve the point_set -- DEPRECATED --
	 * return: multiset of the weighted atoms
	 */
	std::multiset<Weighted_point_3> get_point_set() const { return point_set_; }

private:
	std::vector<double> unit_cell_;
	std::map<Point_3, double> points_;
	std::multiset<Weighted_point_3> point_set_;
	std::vector<std::pair<std::string, Eigen::Vector3d>> atoms_;
	std::vector<double> weights_;
};
