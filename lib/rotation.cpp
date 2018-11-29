// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/rotation.h>
#include <chem/io.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <cassert>
#include <cmath>

Chem::Rotation::Rotation(std::istream& from,
                         const std::string& key,
                         const std::vector<Chem::Element>& at,
                         const Numlib::Mat<double>& x)
    : atms(at), xyz(x), pmom(3), paxis(3, 3), sigma_{1}
{
    using namespace Stdutils;

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "sigma", sigma_, 1);
    }
    rotate_to_principal_axes();
}

void Chem::Rotation::analysis(std::ostream& to) const
{
    if (atms.size() > 1) {
        to << "\nGeometry in principal axes coordinate system:\n";
        Chem::print_geometry(to, atms, xyz);
        Chem::print_center_of_mass(to, center_of_mass());
        Chem::print_principal_moments(to, pmom, paxis);
        Chem::print_rot_constants(to, sigma(), symmetry(), constants());
    }
}

Numlib::Vec<double> Chem::Rotation::constants() const
{
    using namespace Numlib::Constants;

    const double tol = 1.0e-3;
    const double factor = h_bar / (4.0 * pi * giga * m_u * a_0 * a_0 * 1.0e-20);

    Numlib::Vec<double> res(3);

    if (atms.size() > 1) {
        if (std::abs(pmom(0)) < tol) {
            res(0) = factor / pmom(2);
        }
        else {
            res(0) = factor / pmom(0);
            res(1) = factor / pmom(1);
            res(2) = factor / pmom(2);
        }
    }
    return res;
}

std::string Chem::Rotation::symmetry() const
{
    const double tol = 1.0e-3;

    bool ab = std::abs(pmom(0) - pmom(1)) < tol;
    bool bc = std::abs(pmom(1) - pmom(2)) < tol;

    std::string symm = "asymmetric top";
    if (ab && bc) {
        symm = "spherical top";
        if (atms.size() == 1) {
            symm = "atom";
        }
    }
    else if (bc) {
        symm = "prolate symmetric top";
        if (atms.size() == 2) {
            symm = "linear " + symm;
        }
    }
    else if (ab) {
        symm = "oblate symmetric top";
    }
    return symm;
}

Numlib::Vec<double> Chem::Rotation::center_of_mass() const
{
    Numlib::Vec<double> com(3);

    for (Index j = 0; j < xyz.cols(); ++j) {
        double sum = 0.0;
        double tot_mass = 0.0;
        for (Index i = 0; i < xyz.rows(); ++i) {
            sum += atms[i].atomic_mass * xyz(i, j);
            tot_mass += atms[i].atomic_mass;
        }
        com(j) = sum / tot_mass;
    }
    return com;
}

void Chem::Rotation::calc_principal_moments()
{
    // Work on a local copy of the Cartesian coordinates:
    Numlib::Mat<double> xyz_ = xyz;

    // Convert geometry units to bohr:
    xyz_ /= Numlib::Constants::a_0;

    // Move geometry to center of mass:
    auto com = center_of_mass();
    Numlib::translate(xyz_, -com(0), -com(1), -com(2));

    // Compute principal moments:

    paxis = 0.0;
    pmom = 0.0;

    if (atms.size() > 1) {
        for (std::size_t i = 0; i < atms.size(); ++i) {
            double m = atms[i].atomic_mass;
            double x = xyz_(i, 0);
            double y = xyz_(i, 1);
            double z = xyz_(i, 2);
            paxis(0, 0) += m * (y * y + z * z);
            paxis(1, 1) += m * (x * x + z * z);
            paxis(2, 2) += m * (x * x + y * y);
            paxis(0, 1) -= m * x * y;
            paxis(0, 2) -= m * x * z;
            paxis(1, 2) -= m * y * z;
            paxis(1, 0) = paxis(0, 1);
            paxis(2, 0) = paxis(0, 2);
            paxis(2, 1) = paxis(1, 2);
        }
        Numlib::eigs(paxis, pmom);

        // Ensure right-handedness:
        if (Numlib::det(paxis) < 0.0) {
            paxis *= -1.0;
        }
        assert(std::abs(Numlib::det(paxis) - 1.0) < 1.0e-12);
    }
}

