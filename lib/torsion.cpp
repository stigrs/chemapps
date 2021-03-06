// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/torsion.h>
#include <numlib/constants.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cassert>

Chem::Torsion::Torsion(std::istream& from,
                       const std::string& key,
                       const std::vector<Chem::Element>& at,
                       const Numlib::Mat<double>& x,
                       const Numlib::Mat<double>& pa,
                       const Numlib::Vec<double>& pm)
    : atms(at), xyz(x), paxis(pa), pmom(pm)
{
    using namespace Stdutils;

    xyz /= Numlib::Constants::a_0; // convert geometry units to bohr

    perform_analysis = false;

    alpha = Numlib::zeros<Numlib::Mat<double>>(3, 3);

    rot_axis = Numlib::Vec<int>(0);
    rot_top = Numlib::Vec<int>(0);
    sigma_tor = Numlib::Vec<int>(0);

    rmi_tor = Numlib::Vec<double>(0);
    pot_tor = Numlib::Vec<double>(0);
    freq_tor = Numlib::Vec<double>(0);

    x_axis = Numlib::Vec<double>(0);
    y_axis = Numlib::Vec<double>(0);
    z_axis = Numlib::Vec<double>(0);

    top_origo = Numlib::zeros<Numlib::Vec<double>>(3);
    top_com = Numlib::zeros<Numlib::Vec<double>>(3);

    am = 0.0;
    bm = 0.0;
    cm = 0.0;
    um = 0.0;

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "rot_axis", rot_axis, rot_axis);
        get_token_value(from, pos, "rot_top", rot_top, rot_top);
        get_token_value(from, pos, "sigma_tor", sigma_tor, sigma_tor);
        get_token_value(from, pos, "rmi_tor", rmi_tor, rmi_tor);
        get_token_value(from, pos, "pot_tor", pot_tor, pot_tor);
        get_token_value(from, pos, "freq_tor", freq_tor, freq_tor);
    }

    // Calculate reduced moment of inertia if needed:
    if ((rot_top.size() > 0) && rmi_tor.empty()) {
        perform_analysis = true;
        rmi_tor.resize(1);
        rmi_tor(0) = red_moment();
    }

    // Check if data are sensible:
    validate();
}

void Chem::Torsion::analysis(std::ostream& to) const
{
    using namespace Numlib::Constants;

    Stdutils::Format<char> line;
    line.width(24).fill('=');

    Stdutils::Format<double> fix;
    fix.fixed().width(10);

    Stdutils::Format<int> ifix;
    ifix.fixed().width(11);

    Stdutils::Format<double> sci;
    sci.scientific().precision(3);

    if (perform_analysis) {
        to << "\nTorsional Mode Analysis:\n" << line('=') << "\n\n";

        line.width(28).fill('-');
        to << "Atoms defining rotating top:\n"
           << line('-') << '\n'
           << "Center  Atomic  Atomic\n"
           << "Number  Symbol  Mass\n"
           << line('-') << '\n';
        for (Index i = 0; i < rot_top.size(); ++i) {
            to << i + 1 << '\t' << atms[i].atomic_symbol << '\t'
               << fix(atms[i].atomic_mass) << '\n';
        }
        to << line('-') << "\n\n";

        to << "Center " << rot_axis(0) + 1 << " and " << rot_axis(1) + 1
           << " define axis of rotation\n\n";

        fix.width(10).precision(6);
        to << "Center of mass of top (x, y, z): " << fix(top_com(0)) << " "
           << fix(top_com(1)) << " " << fix(top_com(2)) << '\n';

        to << "Origin of coordinates (x, y, z): " << fix(top_origo(0)) << " "
           << fix(top_origo(1)) << " " << fix(top_origo(2)) << "\n\n";

        fix.width(9).precision(3);
        to << "xz product of inertia: " << fix(bm) << " amu bohr^2\n"
           << "yz product of inertia: " << fix(cm) << " amu bohr^2\n"
           << "off-balance factor:    " << fix(um) << " amu bohr^2\n\n";

        fix.width(0).precision(3);
        to << "Moment of inertia of top:  " << fix(am) << " amu bohr^2, "
           << sci(am * au_to_kgm2) << " kg m^2\n";

        to << "Reduced moment of inertia: " << fix(rmi_tor(0))
           << " amu bohr^2, " << sci(rmi_tor(0) * au_to_kgm2) << " kg m^2\n\n";

        const double factor = giga / (c_0 * 100.0);

        to << "Rotational constant: " << fix(constant()(0)) << " GHz, "
           << fix(constant()(0) * factor) << " cm^-1\n";
    }
    else {
        if (!sigma_tor.empty()) {
            line.width(32 + 11 * sigma_tor.size()).fill('-');
            to << "\nTorsional modes:\n" << line('-') << '\n';
            line.width(32).fill(' ');
            to << line(' ');
            for (Index i = 0; i < sigma_tor.size(); ++i) {
                to << "  Minimum " << i + 1;
            }
            line.width(32 + 11 * sigma_tor.size()).fill('-');
            to << '\n' << line('-') << '\n';

            to << "Symmetry number:                ";
            for (auto si : sigma_tor) {
                to << ifix(si);
            }
            to << '\n';

            fix.width(11).precision(3);
            to << "Moment of inertia [amu bohr^2]: ";
            for (auto ri : rmi_tor) {
                to << fix(ri);
            }
            to << '\n';

            to << "Potential energy [cm^-1]:       ";
            for (auto vi : pot_tor) {
                to << fix(vi);
            }
            to << '\n';

            to << "Vibrational frequency [cm^-1]:  ";
            for (auto vi : freq_tor) {
                to << fix(vi);
            }
            to << '\n' << line('-') << '\n';

            to << "Total number of minima:      " << tot_minima() << '\n'
               << "Effective symmetry number:   " << symmetry_number() << '\n'
               << "Effective moment of inertia: " << eff_moment()
               << " amu bohr^2\n";
        }
    }
}

double Chem::Torsion::red_moment()
{
    // Set up axis system for rotating top:

    axis_system();

    // Set up direction cosines matrix:

    direction_cosines();

    // Compute projection of vector from C.O.M. of the molecule to the
    // origin of the coordinates of the top onto the principal axes:

    Numlib::Vec<double> rm(3);
    for (Index i = 0; i < rm.size(); ++i) {
        rm(i) = Numlib::dot(top_origo, paxis.row(i));
    }

    // Calculate moment of inertia of rotating top:

    top_moment_of_inertia();

    // Calculrate reduced moment of inertia (eq. 1):

    Numlib::Vec<double> betam(3);

    for (int i = 0; i < 3; ++i) {
        int im1 = i - 1;
        if (im1 < 0) {
            im1 = 2;
        }
        int ip1 = i + 1;
        if (ip1 > 2) {
            ip1 = 0;
        }
        betam(i) = alpha(2, i) * am - alpha(0, i) * bm - alpha(1, i) * cm +
                   um * (alpha(1, im1) * rm(ip1) - alpha(1, ip1) * rm(im1));
    }
    double tot_mass = 0.0;
    for (std::size_t i = 0; i < atms.size(); ++i) {
        tot_mass += atms[i].atomic_mass;
    }
    double lambdam = 0.0;
    for (int i = 0; i < 3; ++i) {
        lambdam += std::pow(alpha(1, i) * um, 2.0) / tot_mass +
                   std::pow(betam(i), 2.0) / pmom(i);
    }
    return am - lambdam;
}

void Chem::Torsion::validate() const
{
    if (!rot_axis.empty()) {
        if (rot_axis.size() != 2) {
            throw std::runtime_error("bad rot_axis size");
        }
        if (rot_axis(1) > narrow_cast<int>(atms.size())) {
            throw std::runtime_error("bad rot_axis");
        }
    }
    if (!rot_top.empty()) {
        if (rot_top.size() > narrow_cast<int>(atms.size())) {
            throw std::runtime_error("bad rot_top size");
        }
        for (auto ri : rot_top) {
            if (ri > narrow_cast<int>(atms.size())) {
                throw std::runtime_error("bad center in rot_top");
            }
        }
    }
    if (!sigma_tor.empty()) {
        if (Numlib::min(sigma_tor) < 1) {
            throw std::runtime_error("bad sigma_tor");
        }
    }
    if (!rmi_tor.empty()) {
        if (Numlib::min(rmi_tor) <= 0.0) {
            throw std::runtime_error("bad rmi_tor");
        }
    }
    if (!pot_tor.empty()) {
        if (Numlib::min(pot_tor) < 0.0) {
            throw std::runtime_error("bad pot_tor");
        }
    }
    if (!freq_tor.empty()) {
        if (Numlib::min(freq_tor) <= 0.0) {
            throw std::runtime_error("bad freq_tor");
        }
    }
}

void Chem::Torsion::axis_system()
{
    /*
      Definition of axis system for rotating top:

            (x)       z axis is the rotating axis
            COM       x axis is perpendicular to z and goes through the
            /|          center of mass of top
           / |        y axis is perpendicular to x and z axes
        r /  |
         /   |        C1 is the first atom center in rotating axis
        /    |        r is the vector from C1 to center of mass
       /     |_
      /      | |
      C1-------------C2 (z)
    */

    // Calculate center of mass of rotating top and work on a local copy.

    center_of_mass();
    auto com = top_com;

    // Set up z axis - the rotation axis - and its norm:

    z_axis = xyz.row(rot_axis(1)) - xyz.row(rot_axis(0));
    double z_norm = Numlib::norm(z_axis);

    // Find the vector from C1 to C.O.M.:

    auto r_vec = top_com - xyz.row(rot_axis(0));
    double r_norm = Numlib::norm(r_vec);

    // Project r vector onto z axis and find the intersection point:

    const double tol = 1.0e-12;
    double theta = Numlib::dot(r_vec, z_axis) / (r_norm * z_norm);
    if ((std::abs(theta) - 1.0) < tol) { // r and z are parallel
        top_com = xyz.row(rot_top(0));
        r_vec = top_com - xyz.row(rot_axis(0));
        r_norm = Numlib::norm(r_vec);
        theta = Numlib::dot(r_vec, z_axis) / (r_norm * z_norm);
    }
    top_origo = xyz.row(rot_axis(0)) + theta * z_axis * (r_norm / z_norm);

    // Set up x axis:

    x_axis = top_com - top_origo;
    double x_norm = Numlib::norm(x_axis);

    // Check if x and z axes are perpendicular:

    assert(std::abs(Numlib::dot(x_axis, z_axis) / (x_norm * z_norm)) < tol);

    // Generate y axis perpendicular to x and z axes:

    y_axis = Numlib::cross(z_axis, x_axis);
    double y_norm = Numlib::norm(y_axis);

    // Finalize by generating unit vectors:

    x_axis /= x_norm;
    y_axis /= y_norm;
    z_axis /= z_norm;
}

void Chem::Torsion::center_of_mass()
{
    double tot_mass = 0.0;
    for (std::size_t i = 0; i < atms.size(); ++i) {
        tot_mass += atms[i].atomic_mass;
    }
    for (Index j = 0; j < xyz.cols(); ++j) {
        double sum = 0.0;
        for (Index i = 0; i < rot_top.size(); ++i) {
            sum += atms[i].atomic_mass * xyz(i, j);
        }
        top_com(j) = sum / tot_mass;
    }
}

void Chem::Torsion::direction_cosines()
{
    for (int i = 0; i < 3; ++i) {
        alpha(0, i) = Numlib::dot(x_axis, paxis.column(i));
        alpha(1, i) = Numlib::dot(y_axis, paxis.column(i));
        alpha(2, i) = Numlib::dot(z_axis, paxis.column(i));
    }
    if (Numlib::det(alpha) < 0.0) {
        alpha *= -1.0;
    }
    assert(std::abs(Numlib::det(alpha) - 1.0) < 1.0e-12);
}

void Chem::Torsion::top_moment_of_inertia()
{
    // Project coordinates of rotating top onto the (x,y,z) coordinate system:

    Numlib::Mat<double> top_xyz(rot_top.size(), 3);

    for (Index i = 0; i < rot_top.size(); ++i) {
        top_xyz.row(i) = xyz.row(rot_top(i)) - top_origo;
    }
    for (Index i = 0; i < top_xyz.rows(); ++i) {
        double x = Numlib::dot(top_xyz.row(i), x_axis);
        double y = Numlib::dot(top_xyz.row(i), y_axis);
        double z = Numlib::dot(top_xyz.row(i), z_axis);
        top_xyz(i, 0) = x;
        top_xyz(i, 1) = y;
        top_xyz(i, 2) = z;
    }

    // Moment of inertia, products of inertia, and off-balance factor:

    am = 0.0;
    bm = 0.0;
    cm = 0.0;
    um = 0.0;

    for (Index i = 0; i < top_xyz.rows(); ++i) {
        double mass = atms[rot_top(i)].atomic_mass;
        double xi = top_xyz(i, 0);
        double yi = top_xyz(i, 1);
        double zi = top_xyz(i, 2);
        am += mass * (xi * xi + yi * yi);
        bm += mass * xi * zi;
        cm += mass * yi * zi;
        um += mass * xi;
    }
}

double Chem::Torsion::eff_moment() const
{
    double imom_eff = 0.0;
    if (tot_minima() > 0) {
        for (Index i = 0; i < sigma_tor.size(); ++i) {
            imom_eff += sigma_tor(i) * rmi_tor(i); // eq. 7 in C&T (2000)
        }
        imom_eff /= tot_minima();
    }
    return imom_eff;
}

Numlib::Vec<double> Chem::Torsion::constant() const
{
    using namespace Numlib::Constants;

    Numlib::Vec<double> rotc(rmi_tor.size(), 0.0);
    if (!rmi_tor.empty()) {
        for (Index i = 0; i < rotc.size(); ++i) {
            rotc(i) = h_bar / (4.0 * pi * giga * m_u * a_0 * a_0 * 1.0e-20 *
                               rmi_tor(i));
        }
    }
    return rotc;
}

