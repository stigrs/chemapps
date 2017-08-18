/**
   @file torsion.h

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

#ifndef CHEM_TORSION_H
#define CHEM_TORSION_H

#include <chem/molrot.h>
#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <string>

//-----------------------------------------------------------------------------

// Error reporting:

struct Torsion_error : std::runtime_error {
    Torsion_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

/**
   Class for handling torsional modes in molecules using the CT-Cw scheme.

   Algorithm:
   ----------
   The reduced moment of inertia of a symmetrical or unsymmetrical rotating
   top attached to a rigid frame is calculated according to eq 1 in
   Pitzer, K. S. J. Chem. Phys. 1946, vol. 14, pp. 239-243, also known as
   the curvilinear (C) scheme.

   The CT-Cw scheme is described in the following paper: Chuang, Y.-Y.;
   Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.

   Note:
   -----
   Atoms specifying the rotational axis must not be included in the list
   of atoms specifying the rotating top.
*/
class Torsion {
public:
    Torsion(Molrot& rot_) : rot(rot_) {}
    Torsion(std::istream& from, const std::string& key, Molrot& rot_);

    Torsion(const Torsion& tor);

    ~Torsion() {}

    /// Perform torsional mode analysis.
    void analysis(std::ostream& to = std::cout);

    /// Get total number of minima.
    int tot_minima() const;

    /// Calculate effective symmetry number.
    double symmetry_number() const;

    /// Calculate reduced moment of inertia.
    double red_moment_of_inertia();

    /// Calculate effective moment of inertia.
    double eff_moment_of_inertia() const;

    /// Calculate rotational constant for torsional mode.
    arma::vec constant();

private:
    /// Initialize input data.
    void init(std::istream& from, const std::string& key);

    /// Validate input data.
    void validate() const;

    /// Set up axis system for rotating top.
    void axis_system();

    /// Calculate center of mass for rotating top.
    void center_of_mass();

    /// Set up direction cosines matrix.
    void direction_cosines();

    /// Calculate moment of inertia of rotating top.
    void top_moment_of_inertia();

    arma::mat xyz_  = arma::zeros<arma::mat>(0);     ///< local copy of xyz
    arma::mat alpha = arma::zeros<arma::mat>(3, 3);  ///< direction cosines

    arma::uvec rot_axis  = arma::zeros<arma::uvec>(2);  ///< rotational axis
    arma::uvec rot_top   = arma::zeros<arma::uvec>(0);  ///< rotating top moiety
    arma::uvec sigma_tor = arma::zeros<arma::uvec>(0);  ///< symmetry number

    arma::vec rmi_tor  = arma::zeros<arma::vec>(0);  ///< red. moment of inertia
    arma::vec pot_tor  = arma::zeros<arma::vec>(0);  ///< potential coefficients
    arma::vec freq_tor = arma::zeros<arma::vec>(0);  ///< torsional frequencies

    arma::rowvec x_axis = arma::zeros<arma::rowvec>(3);  ///< x axis of rot. top
    arma::rowvec y_axis = arma::zeros<arma::rowvec>(3);  ///< y axis of rot. top
    arma::rowvec z_axis = arma::zeros<arma::rowvec>(3);  ///< z axis of rot. top

    arma::rowvec top_origo = arma::zeros<arma::rowvec>(3);  ///< origo
    arma::rowvec top_com   = arma::zeros<arma::rowvec>(3);  ///< center of mass

    double am;  ///< moment of inertia of rotating top
    double bm;  ///< xz product of inertia
    double cm;  ///< yz product of inertia
    double um;  ///< off-balance factor

    bool perform_torsional_analysis = false;

    Molrot& rot;

    Torsion& operator=(const Torsion& tor);  // no assignments
};

inline Torsion::Torsion(std::istream& from,
                        const std::string& key,
                        Molrot& rot_)
    : rot(rot_)
{
    init(from, key);
}

inline double Torsion::symmetry_number() const
{
    // Eq. 8 in Chuang and Truhlar (2000):
    return tot_minima() / static_cast<double>(sigma_tor.n_elem);
}

#endif /* CHEM_TORSION_H */
