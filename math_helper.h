#pragma once

//--- Boost Includes
#include <boost/math/special_functions.hpp>

//--- Eigen Includes
#include <Eigen/Eigen>

//--- Math Defines
#define PI boost::math::constants::pi<double>()

namespace math {
	/*
	 * Mulitply 3 given values
	 * param1-3: a double value
	 * return: product of the multiplication of the given input
	 */
	inline double mult3(double x, double y, double z) {
		return x * y * z;
	}

	/*
	 * Convert a unit cell from fractional to orthogonal coordinates
	 * param1: vector containing the unit cell lenghts and angles
	 * return: a matrix represeting the converted orthogonal space
	 */
	inline Eigen::Matrix3d fr2cart(std::vector<double> vals) {
		// a, b, c, alpha, beta, gamma in this order from 0 to 5
		Eigen::Matrix3d nm(3, 3);
		
		auto a = vals[0];
		auto b = vals[1];
		auto c = vals[2];
		auto alpha = vals[3];
		auto beta = vals[4];
		auto gamma = vals[5];

		auto cg = boost::math::cos_pi(gamma / 180);
		auto cb = boost::math::cos_pi(beta / 180);
		auto ca = boost::math::cos_pi(alpha / 180);
		auto sg = boost::math::sin_pi(gamma / 180);

		auto v = a * b * c * sqrt(1 - pow(ca, 2) - pow(cb, 2) - pow(cg, 2) + 2 * ca * cb * cg);

		nm(0, 0) = a;
		nm(0, 1) = b * cg;
		nm(0, 2) = c * cb;

		nm(1, 0) = 0;
		nm(1, 1) = b * sg;
		nm(1, 2) = c * (ca - cb * cg) / sg;

		nm(2, 0) = 0;
		nm(2, 1) = 0;
		nm(2, 2) = v / (a * b * sg);

		return nm;
	}
}