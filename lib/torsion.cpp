////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

#include <chem/torsion.h>
#include <srs/datum.h>
#include <srs/utils.h>
#include <algorithm>
#include <cmath>
#include <gsl/gsl>
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
    srs::Format<char> line;
    line.width(24).fill('=');

    srs::Format<double> fix;
    fix.fixed().width(10);

    srs::Format<int> ifix;
    ifix.fixed().width(11);

    srs::Format<double> sci;
    sci.scientific().precision(3);

    if (perform_torsional_analysis) {
        to << "\nTorsional Mode Analysis:\n" << line('=') << "\n\n";

        line.width(28).fill('-');
        to << "Atoms defining rotating top:\n"
           << line('-') << '\n'
           << "Center  Atomic  Atomic\n"
           << "Number  Symbol  Mass\n"
           << line('-') << '\n';
        for (srs::size_t i = 0; i < rot_top.size(); ++i) {
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
            for (srs::size_t i = 0; i < sigma_tor.size(); ++i) {
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
            for (auto pi : pot_tor) {
                to << fix(pi);
            }
            to << '\n';

            to << "Vibrational frequency [cm^-1]:  ";
            for (auto vi : freq_tor) {
                to << fix(vi);
            }
            to << '\n' << line('-') << '\n';

            to << "Total number of minima:      " << tot_minima() << '\n'
               << "Effective symmetry number:   " << symmetry_number() << '\n'
               << "Effective moment of inertia: " << eff_moment_of_inertia()
               << " amu bohr^2\n";
        }
    }
}

double Torsion::eff_moment_of_inertia() const
{
    double imom_eff = 0.0;
    if (tot_minima() > 0) {
        for (srs::size_t i = 0; i < sigma_tor.size(); ++i) {
            imom_eff += sigma_tor(i) * rmi_tor(i);  // eq. 7 in C&T (2000)
        }
        imom_eff /= static_cast<double>(tot_minima());
    }
    return imom_eff;
}

srs::dvector Torsion::constant() const
{
    using namespace datum;

    srs::dvector rotc(rmi_tor.size(), 0.0);
    if (!rmi_tor.empty()) {
        for (srs::size_t i = 0; i < rotc.size(); ++i) {
            rotc(i)
                = h_bar
                  / (4.0 * pi * giga * m_u * a_0 * a_0 * 1.0e-20 * rmi_tor(i));
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

    srs::dvector rm(3);
    for (srs::size_t i = 0; i < rm.size(); ++i) {
        rm(i) = srs::dot(top_origo, rot.paxis.row(i));
    }

    // Calculate moment of inertia of rotating top:

    top_moment_of_inertia();

    // Calculrate reduced moment of inertia (eq. 1):

    srs::dvector betam(3);

    for (int i = 0; i < 3; ++i) {
        int im1 = i - 1;
        if (im1 < 0) {
            im1 = 2;
        }
        int ip1 = i + 1;
        if (ip1 > 2) {
            ip1 = 0;
        }
        betam(i) = alpha(2, i) * am - alpha(0, i) * bm - alpha(1, i) * cm
                   + um * (alpha(1, im1) * rm(ip1) - alpha(1, ip1) * rm(im1));
    }
    double lambdam = 0.0;
    for (int i = 0; i < 3; ++i) {
        lambdam += std::pow(alpha(1, i) * um, 2.0) / rot.tot_mass()
                   + std::pow(betam(i), 2.0) / rot.pmom(i);
    }
    return am - lambdam;
}

void Torsion::init(std::istream& from, const std::string& key)
{
    // Read input data:
    srs::ivector rot_axis_def(2);
    srs::ivector rot_top_def;
    srs::ivector sigma_tor_def;
    srs::dvector rmi_tor_def;
    srs::dvector pot_tor_def;
    srs::dvector freq_tor_def;

    std::map<std::string, srs::Input> input_data;
    input_data["rot_axis"]  = srs::Input(rot_axis, rot_axis_def);
    input_data["rot_top"]   = srs::Input(rot_top, rot_top_def);
    input_data["sigma_tor"] = srs::Input(sigma_tor, sigma_tor_def);
    input_data["rmi_tor"]   = srs::Input(rmi_tor, rmi_tor_def);
    input_data["pot_tor"]   = srs::Input(pot_tor, pot_tor_def);
    input_data["freq_tor"]  = srs::Input(freq_tor, freq_tor_def);

    if (srs::find_section(from, key)) {
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
    for (auto& it : input_data) {
        if (!it.second.is_init()) {
            throw Torsion_error(it.first + " not initialized");
        }
    }

    // Calculate reduced moment of inertia if needed:
    if ((rot_top.size() > 0) && rmi_tor.empty()) {
        perform_torsional_analysis = true;
        rmi_tor.resize(1);
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
    if (rot_axis(1) > gsl::narrow_cast<int>(rot.atoms.size())) {
        throw Torsion_error("bad rot_axis");
    }
    if (rot_top.size() > gsl::narrow_cast<int>(rot.atoms.size())) {
        throw Torsion_error("bad rot_top size");
    }
    for (auto ri : rot_top) {
        if (ri > gsl::narrow_cast<int>(rot.atoms.size())) {
            throw Torsion_error("bad center in rot_top");
        }
    }
    if (!sigma_tor.empty()) {
        if (srs::min(sigma_tor) < 1) {
            throw Torsion_error("bad sigma_tor");
        }
    }
    if (!rmi_tor.empty()) {
        if (srs::min(rmi_tor) <= 0.0) {
            throw Torsion_error("bad rmi_tor");
        }
    }
    if (!pot_tor.empty()) {
        if (srs::min(pot_tor) < 0.0) {
            throw Torsion_error("bad pot_tor");
        }
    }
    if (!freq_tor.empty()) {
        if (srs::min(freq_tor) <= 0.0) {
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
    srs::dvector top_com_ = top_com;

    // Set up z axis - the rotation axis - and its norm:

    z_axis        = xyz_.row(rot_axis(1)) - xyz_.row(rot_axis(0));
    double z_norm = srs::norm(z_axis);

    // Find the vector from C1 to C.O.M.:

    srs::dvector r_vec = top_com_ - xyz_.row(rot_axis(0));
    double r_norm      = srs::norm(r_vec);

    // Project r vector onto z axis and find the intersection point:

    const double tol = 1.0e-12;
    double theta     = srs::dot(r_vec, z_axis) / (r_norm * z_norm);
    if ((std::abs(theta) - 1.0) < tol) {  // r and z are parallel
        top_com_ = xyz_.row(rot_top(0));
        r_vec    = top_com_ - xyz_.row(rot_axis(0));
        r_norm   = srs::norm(r_vec);
        theta    = srs::dot(r_vec, z_axis) / (r_norm * z_norm);
    }
    top_origo = xyz_.row(rot_axis(0)) + theta * z_axis * (r_norm / z_norm);

    // Set up x axis:

    x_axis        = top_com_ - top_origo;
    double x_norm = srs::norm(x_axis);

    // Check if x and z axes are perpendicular:

    double xz_angle = srs::dot(x_axis, z_axis) / (x_norm * z_norm);
    Expects(std::abs(xz_angle) < tol);

    // Generate y axis perpendicular to x and z axes:

    y_axis        = srs::cross(z_axis, x_axis);
    double y_norm = srs::norm(y_axis);

    // Finalize by generating unit vectors:

    x_axis /= x_norm;
    y_axis /= y_norm;
    z_axis /= z_norm;
}

void Torsion::center_of_mass()
{
    std::cout << top_com.size() << std::endl;
    for (srs::size_t j = 0; j < xyz_.cols(); ++j) {
        double sum = 0.0;
        for (srs::size_t i = 0; i < rot_top.size(); ++i) {
            sum += rot.atoms[i].atomic_mass * xyz_(i, j);
        }
        top_com(j) = sum / rot.tot_mass();
    }
}

void Torsion::direction_cosines()
{
    for (int i = 0; i < 3; ++i) {
        alpha(0, i) = srs::dot(x_axis, rot.paxis.column(i));
        alpha(1, i) = srs::dot(y_axis, rot.paxis.column(i));
        alpha(2, i) = srs::dot(z_axis, rot.paxis.column(i));
    }
    Ensures(std::abs(srs::det(alpha) - 1.0) < 1.0e-12);
}

void Torsion::top_moment_of_inertia()
{
    // Project coordinates of rotating top onto the (x,y,z) coordinate system:

    srs::dmatrix top_xyz(rot_top.size(), 3);
    for (srs::size_t i = 0; i < rot_top.size(); ++i) {
        top_xyz.row(i) = xyz_.row(rot_top(i)) - top_origo;
    }
    for (srs::size_t i = 0; i < top_xyz.rows(); ++i) {
        double x      = srs::dot(top_xyz.row(i), x_axis);
        double y      = srs::dot(top_xyz.row(i), y_axis);
        double z      = srs::dot(top_xyz.row(i), z_axis);
        top_xyz(i, 0) = x;
        top_xyz(i, 1) = y;
        top_xyz(i, 2) = z;
    }

    // Moment of inertia, products of inertia, and off-balance factor:

    am = 0.0;
    bm = 0.0;
    cm = 0.0;
    um = 0.0;

    for (int i = 0; i < top_xyz.rows(); ++i) {
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
