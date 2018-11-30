// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_TORSION_H
#define CHEM_TORSION_H

#include <chem/element.h>
#include <numlib/traits.h>
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

// Class for handling torsional modes in molecules using the CT-Cw scheme.
//
// Algorithm:
//   The reduced moment of inertia of a symmetrical or unsymmetrical rotating
//   top attached to a rigid frame is calculated according to eq 1 in
//   Pitzer, K. S. J. Chem. Phys. 1946, vol. 14, pp. 239-243, also known as
//   the curvilinear (C) scheme.
//
//   The CT-Cw scheme is described in the following paper: Chuang, Y.-Y.;
//   Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.
//
// Note:
//   Atoms specifying the rotational axis must not be included in the list
//   of atoms specifying the rotating top.
//
class Torsion {
public:
    Torsion() = default;

    Torsion(std::istream& from,
            const std::string& key,
            const std::vector<Element>& at,
            const Numlib::Mat<double>& x,
            const Numlib::Mat<double>& pa,
            const Numlib::Vec<double>& pm);

    // Copy semantics:
    Torsion(const Torsion&) = default;
    Torsion& operator=(const Torsion&) = default;

    // Move semantics:
    Torsion(Torsion&&) = default;
    Torsion& operator=(Torsion&&) = default;

    ~Torsion() = default;

    // Perform torsional mode analysis.
    void analysis(std::ostream& to) const;

    // Set new Cartesian coordinates.
    void set(const Numlib::Mat<double>& x,
             const Numlib::Mat<double>& pa,
             const Numlib::Vec<double>& pm);

    // Get total number of minima (eq 1 in C&T, 2000).
    int tot_minima() const { return Numlib::sum(sigma_tor); }

    // Calculate effective symmetry number.
    double symmetry_number() const;

    // Calculate reduced moment of inertia.
    double red_moment();

    // Calculate effective moment of inertia.
    double eff_moment() const;

    // Calculate rotational constant for torsional mode.
    Numlib::Vec<double> constant() const;

    // Return potential coefficients.
    const auto& pot_coeff() const { return pot_tor; }

    // Return torsional frequencies.
    const auto& frequencies() const { return freq_tor; }

private:
    // Validate input data.
    void validate() const;

    // Set up axis system for rotating top.
    void axis_system();

    // Calculate center of mass for rotating top.
    void center_of_mass();

    // Set up direction cosines matrix.
    void direction_cosines();

    // Calculate moment of inertia of rotating top.
    void top_moment_of_inertia();

    std::vector<Element> atms;

    Numlib::Mat<double> xyz;
    Numlib::Mat<double> paxis;
    Numlib::Vec<double> pmom;

    Numlib::Mat<double> alpha; // direction cosines

    Numlib::Vec<int> rot_axis;  // rotational axis
    Numlib::Vec<int> rot_top;   // rotating top moiety
    Numlib::Vec<int> sigma_tor; // symmetry number

    Numlib::Vec<double> rmi_tor;  // red. moment of inertia
    Numlib::Vec<double> pot_tor;  // potential coefficients
    Numlib::Vec<double> freq_tor; // torsional frequencies

    Numlib::Vec<double> x_axis; // x axis of rotating top
    Numlib::Vec<double> y_axis; // y axis of rotating top
    Numlib::Vec<double> z_axis; // z axis of rotating top

    Numlib::Vec<double> top_origo; // origo of rotating top
    Numlib::Vec<double> top_com;   // center of mass of rotating top

    double am; // moment of inertia of rotating top
    double bm; // xz product of inertia
    double cm; // yz product of inertia
    double um; // off-balance factor

    bool perform_analysis;
};

inline void Torsion::set(const Numlib::Mat<double>& x,
                         const Numlib::Mat<double>& pa,
                         const Numlib::Vec<double>& pm)
{
    Assert::dynamic(Numlib::same_extents(xyz, x),
                    "bad size of Cartesian coordinates");
    Assert::dynamic(Numlib::same_extents(paxis, pa),
                    "bad size of principal axes");
    Assert::dynamic(Numlib::same_extents(pmom, pm),
                    "bad size of principal moments");
    xyz = x;
    paxis = pa;
    pmom = pm;
}

inline double Torsion::symmetry_number() const
{
    // Eq. 8 in Chuang and Truhlar (2000):
    return tot_minima() / narrow_cast<double>(sigma_tor.size());
}

} // namespace Chem

#endif // CHEM_TORSION_H

