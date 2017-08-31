/**
   @file math.h

   This file is part of ChemApps - A C++ Chemistry Toolkit

   Copyright (C) 2016-2017  Stig Rune Sellevag

   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CHEM_MATH_H
#define CHEM_MATH_H

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4100)  // caused by armadillo
#endif                           // _MSC_VER

#include <chem/datum.h>
#include <armadillo>
#include <cmath>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER

namespace chem {

// Convert radians to degrees.
inline double radtodeg(double rad) { return rad * 180.0 / datum::pi; }

// Convert degrees to radians.
inline double degtorad(double deg) { return deg * datum::pi / 180.0; }

// Compute distance between two points.
double distance(const arma::vec& a, const arma::vec& b);

// Compute angle in degrees between three points.
double angle(const arma::vec& a, const arma::vec& b, const arma::vec& c);

// Compute dihedral angle in degrees given four points.
double dihedral(const arma::vec& a,
                const arma::vec& b,
                const arma::vec& c,
                const arma::vec& d);

// Compute the pair-wise distances between observations in n-dim. space.
void pdist_matrix(arma::mat& dm, const arma::mat& mat);

// Perform translation.
void translate(arma::mat& xyz, double dx, double dy, double dz);

// Perform rotation given a rotation matrix.
void rotate(arma::mat& xyz, const arma::mat33& rotm);

// Compute abscissas and weights of Gauss-Legendre n-point quadrature formula.
void gaussleg(
    int n, arma::vec& x, arma::vec& w, double a = -1.0, double b = 1.0);

// Floating point comparison.
bool approx_equal(double a, double b, double epsilon = 1.0e-12);

// Check if integer is even.
bool is_even(int n) { return n % 2 ? false : true; }

// Check if integer is odd.
bool is_odd(int n) { return n % 2 ? true : false; }

}  // namespace chem

inline double chem::distance(const arma::vec& a, const arma::vec& b)
{
    return arma::norm(b - a);
}

inline double chem::angle(const arma::vec& a,
                          const arma::vec& b,
                          const arma::vec& c)
{
    arma::vec ab = arma::normalise(a - b);
    arma::vec bc = arma::normalise(c - b);
    return radtodeg(std::acos(arma::dot(ab, bc)));
}

inline bool chem::approx_equal(double a, double b, double epsilon)
{
    using namespace std;
    return abs(a - b) <= ((abs(a) < abs(b) ? abs(b) : abs(a)) * epsilon);
}

#endif /* CHEM_MATH_H */
