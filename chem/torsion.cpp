//////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
///////////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4100)  // caused by armadillo
#endif                           // _MSC_VER

#include <chem/datum.h>
#include <chem/input.h>
#include <chem/torsion.h>
#include <chem/utils.h>
#include <cmath>
#include <iostream>
#include <map>

Torsion::Torsion(const Torsion& tor) : rot(tor.rot)
{
    xyz_  = tor.xyz_;
    alpha = tor.alpha;

    rot_axis  = tor.rot_axis;
    rot_top   = tor.rot_top;
    sigma_tor = tor.sigma_tor;

    rmi_tor  = tor.rmi_tor;
    pot_tor  = tor.pot_tor;
    freq_tor = tor.freq_tor;

    x_axis = tor.x_axis;
    y_axis = tor.y_axis;
    z_axis = tor.z_axis;

    top_origo = tor.top_origo;
    top_com   = tor.top_com;

    am = tor.am;
    bm = tor.bm;
    cm = tor.cm;
    um = tor.um;

    perform_torsional_analysis = tor.perform_torsional_analysis;
}

void Torsion::analysis(std::ostream& to)
{
    chem::Format<char> line;
    line.width(24).fill('=');

    chem::Format<double> fix;
    fix.fixed().width(10);

    chem::Format<int> ifix;
    ifix.fixed().width(11);

    chem::Format<double> sci;
    sci.scientific().precision(3);

    if (perform_torsional_analysis) {
        to << "\nTorsional Mode Analysis:\n" << line('=') << "\n\n";

        line.width(28).fill('-');
        to << "Atoms defining rotating top:\n"
           << line('-') << '\n'
           << "Center  Atomic  Atomic\n"
           << "Number  Symbol  Mass\n"
           << line('-') << '\n';
        for (arma::uword i = 0; i < rot_top.size(); ++i) {
            to << i + 1 << '\t' << rot.atoms[i].atomic_symbol << '\t'
               << fix(rot.atoms[i].atomic_mass) << '\n';
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
           << sci(am * datum::au_to_kgm2) << " kg m^2\n";

        to << "Reduced moment of inertia: " << fix(rmi_tor(0))
           << " amu bohr^2, " << sci(rmi_tor(0) * datum::au_to_kgm2)
           << " kg m^2\n\n";

        const double factor = datum::giga / (datum::c_0 * 100.0);

        to << "Rotational constant: " << fix(constant()(0)) << " GHz, "
           << fix(constant()(0) * factor) << " cm^-1\n";
    }
    else {
        if (!sigma_tor.empty()) {
            line.width(32 + 11 * sigma_tor.size()).fill('-');
            to << "\nTorsional modes:\n" << line('-') << '\n';
            line.width(32).fill(' ');
            to << line(' ');
            for (arma::uword i = 0; i < sigma_tor.size(); ++i) {
                to << "  Minimum " << i + 1;
            }
            line.width(32 + 11 * sigma_tor.size()).fill('-');
            to << '\n' << line('-') << '\n';

            to << "Symmetry number:                ";
            for (arma::uword i = 0; i < sigma_tor.size(); ++i) {
                to << ifix(sigma_tor(i));
            }
            to << '\n';

            fix.width(11).precision(3);
            to << "Moment of inertia [amu bohr^2]: ";
            for (arma::uword i = 0; i < rmi_tor.size(); ++i) {
                to << fix(rmi_tor(i));
            }
            to << '\n';

            to << "Potential energy [cm^-1]:       ";
            for (arma::uword i = 0; i < pot_tor.size(); ++i) {
                to << fix(pot_tor(i));
            }
            to << '\n';

            to << "Vibrational frequency [cm^-1]:  ";
            for (arma::uword i = 0; i < freq_tor.size(); ++i) {
                to << fix(freq_tor(i));
            }
            to << '\n' << line('-') << '\n';

            to << "Total number of minima:      " << tot_minima() << '\n'
               << "Effective symmetry number:   " << symmetry_number() << '\n'
               << "Effective moment of inertia: " << eff_moment_of_inertia()
               << " amu bohr^2\n";
        }
    }
}

int Torsion::tot_minima() const
{
    int nminima = 0;
    for (arma::uword i = 0; i < sigma_tor.size(); ++i) {
        nminima += sigma_tor(i);  // eq. 1 in C&T (2000)
    }
    return nminima;
}

double Torsion::eff_moment_of_inertia() const
{
    double imom_eff = 0.0;
    if (tot_minima() > 0) {
        for (arma::uword i = 0; i < sigma_tor.size(); ++i) {
            imom_eff += sigma_tor(i) * rmi_tor(i);  // eq. 7 in C&T (2000)
        }
        imom_eff /= static_cast<double>(tot_minima());
    }
    return imom_eff;
}

arma::vec Torsion::constant()
{
    using namespace datum;

    arma::vec rotc = arma::zeros<arma::vec>(rmi_tor.size());
    if (!rmi_tor.empty()) {
        for (arma::uword i = 0; i < rotc.size(); ++i) {
            rotc(i) = h_bar / (4.0 * pi * giga * m_u * a_0 * a_0 * 1.0e-20 *
                               rmi_tor(i));
        }
    }
    return rotc;
}

double Torsion::red_moment_of_inertia()
{
    // Rotate molecule to principal axes and compute principal moments:

    rot.rotate_to_principal_axes();
    rot.principal_moments();

    // Work on a local copy of the XYZ and convert coordinates to bohr:

    xyz_ = rot.xyz;
    xyz_ /= datum::a_0;

    // Set up axis system for rotating top:

    axis_system();

    // Set up direction cosines matrix:

    direction_cosines();

    // Compute projection of vector from C.O.M. of the molecule to the
    // origin of the coordinates of the top onto the principal axes:

    arma::rowvec rm = arma::zeros<arma::rowvec>(3);
    for (arma::uword i = 0; i < rm.size(); ++i) {
        rm(i) = arma::dot(top_origo, rot.paxis.row(i));
    }

    // Calculate moment of inertia of rotating top:

    top_moment_of_inertia();

    // Calculrate reduced moment of inertia (eq. 1):

    arma::rowvec betam = arma::zeros<arma::rowvec>(3);

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
    double lambdam = 0.0;
    for (int i = 0; i < 3; ++i) {
        lambdam += std::pow(alpha(1, i) * um, 2.0) / rot.tot_mass() +
                   std::pow(betam(i), 2.0) / rot.pmom(i);
    }
    return am - lambdam;
}

void Torsion::init(std::istream& from, const std::string& key)
{
    // Read input data:
    arma::uvec rot_axis_def  = arma::zeros<arma::uvec>(2);
    arma::uvec rot_top_def   = arma::zeros<arma::uvec>(0);
    arma::uvec sigma_tor_def = arma::zeros<arma::uvec>(0);
    arma::vec rmi_tor_def    = arma::zeros<arma::vec>(0);
    arma::vec pot_tor_def    = arma::zeros<arma::vec>(0);
    arma::vec freq_tor_def   = arma::zeros<arma::vec>(0);

    std::map<std::string, Input> input_data;
    input_data["rot_axis"]  = Input(rot_axis, rot_axis_def);
    input_data["rot_top"]   = Input(rot_top, rot_top_def);
    input_data["sigma_tor"] = Input(sigma_tor, sigma_tor_def);
    input_data["rmi_tor"]   = Input(rmi_tor, rmi_tor_def);
    input_data["pot_tor"]   = Input(pot_tor, pot_tor_def);
    input_data["freq_tor"]  = Input(freq_tor, freq_tor_def);

    if (chem::find_section(from, key)) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }
    else {
        throw Torsion_error("cannot find " + key + " section");
    }

    // Check if initialized:
    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Torsion_error(it->first + " not initialized");
        }
    }

    // Calculate reduced moment of inertia if needed:
    if ((rot_top.size() > 0) && rmi_tor.empty()) {
        perform_torsional_analysis = true;
        rmi_tor.set_size(1);
        rmi_tor(0) = red_moment_of_inertia();
    }

    // Check if data are sensible:
    validate();
}

void Torsion::validate() const
{
    if (rot_axis.size() != 2) {
        throw Torsion_error("bad rot_axis size");
    }
    if (rot_axis(1) > rot.atoms.size()) {
        throw Torsion_error("bad rot_axis");
    }
    if (rot_top.size() > rot.atoms.size()) {
        throw Torsion_error("bad rot_top size");
    }
    for (arma::uword i = 0; i < rot_top.size(); ++i) {
        if (rot_top(i) > rot.atoms.size()) {
            throw Torsion_error("bad center in rot_top");
        }
    }
    if (sigma_tor.size() > 0) {
        if (arma::any(sigma_tor < 1)) {
            throw Torsion_error("bad sigma_tor");
        }
    }
    if (rmi_tor.size() > 0) {
        if (arma::any(rmi_tor <= 0.0)) {
            throw Torsion_error("bad rmi_tor");
        }
    }
    if (pot_tor.size() > 0) {
        if (arma::any(pot_tor < 0.0)) {
            throw Torsion_error("bad pot_tor");
        }
    }
    if (freq_tor.size() > 0) {
        if (arma::any(freq_tor <= 0.0)) {
            throw Torsion_error("bad freq_tor");
        }
    }
}

void Torsion::axis_system()
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
    arma::rowvec top_com_(top_com);

    // Set up z axis - the rotation axis - and its norm:

    z_axis        = xyz_.row(rot_axis(1)) - xyz_.row(rot_axis(0));
    double z_norm = arma::norm(z_axis);

    // Find the vector from C1 to C.O.M.:

    arma::rowvec r_vec = top_com - xyz_.row(rot_axis(0));
    double r_norm      = arma::norm(r_vec);

    // Project r vector onto z axis and find the intersection point:

    const double tol = 1.0e-12;
    double theta     = arma::dot(r_vec, z_axis) / (r_norm * z_norm);
    if ((std::abs(theta) - 1.0) < tol) {  // r and z are parallel
        top_com = xyz_.row(rot_top(0));
        r_vec   = top_com - xyz_.row(rot_axis(0));
        r_norm  = arma::norm(r_vec);
        theta   = arma::dot(r_vec, z_axis) / (r_norm * z_norm);
    }
    top_origo = xyz_.row(rot_axis(0)) + theta * z_axis * r_norm / z_norm;

    // Set up x axis:

    x_axis        = top_com - top_origo;
    double x_norm = arma::norm(x_axis);

    // Check if x and z axes are perpendicular:

    double xz_angle = arma::dot(x_axis, z_axis) / (x_norm * z_norm);
    chem::Assert(std::abs(xz_angle) < tol,
                 Torsion_error("x and z axes are not parallel"));

    // Generate y axis perpendicular to x and z axes:

    y_axis        = arma::cross(z_axis, x_axis);
    double y_norm = arma::norm(y_axis);

    // Finalize by generating unit vectors:

    x_axis /= x_norm;
    y_axis /= y_norm;
    z_axis /= z_norm;
}

void Torsion::center_of_mass()
{
    for (arma::uword j = 0; j < xyz_.n_cols; ++j) {
        double sum = 0.0;
        for (arma::uword i = 0; i < rot_top.size(); ++i) {
            sum += rot.atoms[i].atomic_mass * xyz_(i, j);
        }
        top_com(j) = sum / rot.tot_mass();
    }
}

void Torsion::direction_cosines()
{
    for (arma::uword i = 0; i < 3; ++i) {
        alpha(0, i) = arma::dot(x_axis, rot.paxis.col(i));
        alpha(1, i) = arma::dot(y_axis, rot.paxis.col(i));
        alpha(2, i) = arma::dot(z_axis, rot.paxis.col(i));
    }
    chem::Assert(std::abs(arma::det(alpha) - 1.0) < 1.0e-12,  // should be +1
                 Torsion_error("bad direction cosines matrix"));
}

void Torsion::top_moment_of_inertia()
{
    // Project coordinates of rotating top onto the (x,y,z) coordinate system:

    arma::mat top_xyz(rot_top.size(), 3);
    for (arma::uword i = 0; i < rot_top.size(); ++i) {
        top_xyz.row(i) = xyz_.row(rot_top(i)) - top_origo;
    }
    for (arma::uword i = 0; i < top_xyz.n_rows; ++i) {
        double x = arma::dot(top_xyz.row(i), x_axis);
        double y = arma::dot(top_xyz.row(i), y_axis);
        double z = arma::dot(top_xyz.row(i), z_axis);
        top_xyz(i, 0) = x;
        top_xyz(i, 1) = y;
        top_xyz(i, 2) = z;
    }

    // Moment of inertia, products of inertia, and off-balance factor:

    am = 0.0;
    bm = 0.0;
    cm = 0.0;
    um = 0.0;

    for (arma::uword i = 0; i < top_xyz.n_rows; ++i) {
        double mass = rot.atoms[rot_top(i)].atomic_mass;
        double xi   = top_xyz(i, 0);
        double yi   = top_xyz(i, 1);
        double zi   = top_xyz(i, 2);
        am += mass * (xi * xi + yi * yi);
        bm += mass * xi * zi;
        cm += mass * yi * zi;
        um += mass * xi;
    }
}
