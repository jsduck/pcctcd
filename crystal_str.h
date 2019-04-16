#pragma once

//-- Standard Includes
#include <vector>

//-- Local Includes
#include <Eigen/Eigen>

#include "cgal_helper.h"
#include <unordered_map>

class Crystal
{
public:
	struct d3
	{
		d3(double nx, double ny, double nz) {
			x = nx;
			y = ny;
			z = nz;
		}
		double x, y, z;
	};

public:
	Crystal();
	
	void clear() {
		unit_cell_.clear();
		atoms_.clear();
	}

	void append_cell(std::vector<double> v);
	void reduce_cell(double epsilon, int loops = 50);

	std::vector<double>& get_cell();

	void append_atom(std::string label, Eigen::Vector3d coords);

	std::pair<std::string, Eigen::Vector3d>& get_atom(int index);
	std::vector<std::pair<std::string, Eigen::Vector3d>>& get_atoms();

	void orthogonal();

	std::vector<d3> get_points();
	std::map<Point_3, double> get_orthogonal_points() const;

	std::vector<double> get_vdW_weights() const;

	void output_off(std::string filename);
	void output_xyz(std::string filename);

private:
	std::vector<double> unit_cell_;
	std::map<Point_3, double> points_;
	std::vector<std::pair<std::string, Eigen::Vector3d>> atoms_;

	
};
