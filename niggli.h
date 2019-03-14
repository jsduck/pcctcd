#pragma once

//--- Boost includes
#include <boost/units/cmath.hpp>
#include <boost/math/special_functions.hpp>

struct niggli
{
	static bool reduce(std::vector<double> &vals, double epsilon, int max_loops = 50) {
		using namespace boost::math;

		// variable for breaking out of the reduction loop
		int br = 0;

		const auto a = vals[0];
		const auto b = vals[1];
		const auto c = vals[2];
		const auto alpha = vals[3];
		const auto beta = vals[4];
		const auto gamma = vals[5];

		// printf("INPUT:\n a: %.4f b: %.4f c: %.4f \n alpha: %.4f beta: %.4f gamma: %.4f\n", a, b, c, alpha, beta, gamma);


		//--- Setting up initial values of variables
		double A = pow(a, 2);
		double B = pow(b, 2);
		double C = pow(c, 2);

		double xi = 2 * b * c * cos_pi(alpha / 180.f);
		double eta = 2 * a * c * cos_pi(beta / 180.f);
		double zeta = 2 * a * b * cos_pi(gamma / 180.f);

		// TEST VALUES
		// A = 9; B = 27; C = 4;
		// xi = -5; eta = -4; zeta = -22;

		// Implementation of the "Unified Algorithm for Determining the Reduced (Niggli) Cell" by I. KRIVY and B. GRUBER - 1975
		// STEP 1
		step1:
		if (br > max_loops) {
			//printf("NIGGLI CELL NOT REDUCED\n");
			return false;
		}
		br++;
		if (A > B + epsilon || (!(abs(A - B) > epsilon) && abs(xi) > abs(eta) + epsilon)) {
			std::swap(A, B);
			std::swap(xi, eta);
		}
		// STEP 2
		if (B > C + epsilon || (!(abs(B - C) > epsilon) && abs(eta) > abs(zeta) + epsilon)) {
			std::swap(B, C);
			std::swap(eta, zeta);

			goto step1;
		}
		// STEP 3
		if (xi * eta * zeta > 0) {
			xi = abs(xi);
			eta = abs(eta);
			zeta = abs(zeta);

		}
		// STEP 4
		if (xi * eta * zeta <= 0) {
			xi = -abs(xi);
			eta = -abs(eta);
			zeta = -abs(zeta);

		}
		// STEP 5
		if (abs(xi) > B + epsilon || (!(abs(B - xi) > epsilon) && 2 * eta < zeta - epsilon) || (!(abs(B + xi) > epsilon) && zeta < -epsilon)) {
			C = B + C - xi * sign(xi);
			eta = eta - zeta * sign(xi);
			xi = xi - 2 * B * sign(xi);

			goto step1;
		}
		// STEP 6
		if (abs(eta) > A + epsilon || (!(abs(A - eta) > epsilon) && 2 * xi < zeta - epsilon) || (!(abs(A + eta) > epsilon) && zeta < -epsilon)) {
			C = A + C - eta * sign(eta);
			xi = xi - zeta * sign(eta);
			eta = eta - 2 * A * sign(eta);

			goto step1;
		}
		// STEP 7
		if (abs(zeta) > A + epsilon || (!(abs(A - zeta) > epsilon) && 2 * xi < eta - epsilon) || (!(abs(A + zeta) > epsilon) && eta < -epsilon)) {
			B = A + B - zeta * sign(zeta);
			xi = xi - eta * sign(zeta);
			zeta = zeta - 2 * A * sign(zeta);

			goto step1;
		}
		// STEP 8
		if (xi + eta + zeta + A + B < -epsilon || (!(abs(xi + eta + zeta + A + B) > epsilon) && 2 * (A + eta) + zeta > epsilon)) {
			C = A + B + C + xi + eta + zeta;
			xi = 2 * B + xi + zeta;
			eta = 2 * A + eta + zeta;

			goto step1;
		}

		// printf("NIGGLI CELL REDUCED\n");

		// Normalise the results to describe the Niggli cell in the usual way
		A = sqrt(A);
		B = sqrt(B);
		C = sqrt(C);

		const auto PI = boost::math::constants::pi<double>();

		xi = acos(xi / (B * C) / 2) / PI * 180;
		eta = acos(eta / (A * C) / 2) / PI * 180;
		zeta = acos(zeta / (A * B) / 2) / PI * 180;

		// printf("NORM:\n a': %.4f b': %.4f c': %.4f \n alpha': %.4f beta': %.4f gamma': %.4f\n", A, B, C, xi, eta, zeta);

		vals[0] = A;
		vals[1] = B;
		vals[2] = C;
		vals[3] = xi;
		vals[4] = eta;
		vals[5] = zeta;

		return true;
	}
};

