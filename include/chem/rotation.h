// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_ROTATION_H
#define CHEM_ROTATION_H

#include <chem/element.h>
#include <numlib/math.h>
#include <numlib/matrix.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

// Class for handling molecular rotations.
//
class Rotation {
public:
    Rotation() = default;

    Rotation(const std::vector<Element>& at,
             const Numlib::Mat<double>& x,
             int sig = 1);

    Rotation(std::istream& from,
             const std::string& key,
             const std::vector<Element>& at,
             const Numlib::Mat<double>& x);

    // Copy semantics:
    Rotation(const Rotation&) = default;
    Rotation& operator=(const Rotation&) = default;

    // Move semantics:
    Rotation(Rotation&&) = default;
    Rotation& operator=(Rotation&&) = default;

    ~Rotation() = default;

    // Perform rotational analysis.
    void analysis(std::ostream& to = std::cout) const;

    // Get rotational symmetry number.
    auto sigma() const { return sigma_; }

    // Compute rotational constants.
    Numlib::Vec<double> constants() const;

    // Compute rotational symmetry.
    std::string symmetry() const;

    // Get principal moments.
    auto principal_moments() const { return pmom; }

    // Get principal axes.
    auto principal_axes() const { return paxis; }

    // Get Cartesian coordinates in principal axes coordinate system.
    const auto& get_xyz_paxis() const { return xyz; }

private:
    // Move geometry to center of mass.
    void move_to_com();

    // Rotate to principal axes.
    //
    // Note: This coordinate system is not the same as Gaussian's
    // standard orientation.
    //
    void rotate_to_principal_axes();

    // Compute center of mass coordinates.
    Numlib::Vec<double> center_of_mass() const;

    // Compute principal moments of inertia.
    void calc_principal_moments();

    std::vector<Element> atms;
    Numlib::Mat<double> xyz;

    Numlib::Vec<double> pmom;
    Numlib::Mat<double> paxis;

    int sigma_;
};

inline Rotation::Rotation(const std::vector<Element>& at,
                          const Numlib::Mat<double>& x,
                          int sig)
    : atms(at), xyz(x), pmom(3), paxis(3, 3), sigma_{sig}
{
    rotate_to_principal_axes();
}

inline void Rotation::move_to_com()
{
    auto com = center_of_mass();
    Numlib::translate(xyz, -com(0), -com(1), -com(2));
}

inline void Rotation::rotate_to_principal_axes()
{
    if (!atms.empty()) {
        move_to_com();
        calc_principal_moments();
    }
}

} // namespace Chem

#endif // CHEM_ROTATION_H

