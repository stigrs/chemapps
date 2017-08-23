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

#include <chem/datum.h>
#include <chem/input.h>
#include <chem/molecule_io.h>
#include <chem/molrot.h>
#include <chem/utils.h>
#include <cmath>
#include <map>

Molrot::Molrot(const Molrot& rot)
    : atoms(rot.atoms), xyz(rot.xyz), pmom(rot.pmom), paxis(rot.paxis)
{
    sigma   = rot.sigma;
    aligned = rot.aligned;
}

void Molrot::analysis(std::ostream& to)
{
    if (!atoms.empty()) {
        rotate_to_principal_axes();
        to << "\nGeometry in principal axes coordinate system:\n";
        chem::print_geometry(to, atoms, xyz);
        print_center_of_mass(to);
        print_principal_moments(to);
        print_constants(to);
    }
}

double Molrot::tot_mass() const
{
    double totmass = 0.0;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        totmass += atoms[i].atomic_mass;
    }
    return totmass;
}

arma::vec3 Molrot::constants()
{
    using namespace datum;

    if (!aligned) {
        rotate_to_principal_axes();
    }
    const double tol    = 1.0e-3;
    const double factor = h_bar / (4.0 * pi * giga * m_u * a_0 * a_0 * 1.0e-20);

    arma::vec3 rotc;
    rotc.zeros();

    if (atoms.size() > 1) {
        if (std::abs(pmom(0)) < tol) {
            rotc(0) = factor / pmom(2);
        }
        else {
            rotc(0) = factor / pmom(0);
            rotc(1) = factor / pmom(1);
            rotc(2) = factor / pmom(2);
        }
    }
    return rotc;
}

std::string Molrot::symmetry()
{
    if (!aligned) {
        rotate_to_principal_axes();
    }
    const double tol = 1.0e-3;

    bool ab = std::abs(pmom(0) - pmom(1)) < tol;
    bool bc = std::abs(pmom(1) - pmom(2)) < tol;

    std::string symm = "asymmetric top";
    if (ab && bc) {
        symm = "spherical top";
        if (atoms.size() == 1) {
            symm = "atom";
        }
    }
    else if (bc) {
        symm = "prolate symmetric top";
        if (atoms.size() == 2) {
            symm = "linear " + symm;
        }
    }
    else if (ab) {
        symm = "oblate symmetric top";
    }
    return symm;
}

void Molrot::init(std::istream& from, const std::string& key)
{
    // Read input data:

    std::map<std::string, Input> input_data;
    input_data["sigma"] = Input(sigma, 1.0);

    bool found = chem::find_section(from, key);
    if (found) {
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
        throw Molrot_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Molrot_error(it->first + " not initialized");
        }
    }
}

arma::vec3 Molrot::center_of_mass() const
{
    arma::vec3 com;

    for (arma::uword j = 0; j < xyz.n_cols; ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            sum += atoms[i].atomic_mass * (xyz)(i, j);
        }
        com(j) = sum / tot_mass();
    }
    return com;
}

void Molrot::principal_moments()
{
    // Work on a local copy of the cartesian coordinates:
    arma::mat xyz_(xyz);

    // Convert geometry to bohr:
    xyz_ /= datum::a_0;

    // Move geometry to center of mass:
    arma::vec3 com = center_of_mass();
    chem::translate(xyz_, -com(0), -com(1), -com(2));

    // Compute principal moments:
    paxis.zeros();
    pmom.zeros();

    if (atoms.size() > 1) {
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            double m = atoms[i].atomic_mass;
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
        arma::eig_sym(pmom, paxis, paxis);

        // Check if principal axes form a proper rotation:
        if (arma::det(paxis) < 0.0) {
            paxis *= -1.0;
        }
        chem::Assert((std::abs(arma::det(paxis)) - 1.0) < 1.0e-12,
                     Molrot_error("bad determinant of paxis"));
    }
}

void Molrot::rotate_to_principal_axes()
{
    if (!atoms.empty()) {
        move_to_com();
        principal_moments();
        chem::rotate(xyz, paxis.t());
        principal_moments();
        aligned = true;
    }
}

void Molrot::print_center_of_mass(std::ostream& to) const
{
    arma::vec3 com = center_of_mass();
    chem::Format<double> fix;
    fix.fixed().width(8).precision(4);

    if (!atoms.empty()) {
        to << "Center of mass (X, Y, Z):  " << fix(com(0)) << ", "
           << fix(com(1)) << ", " << fix(com(2)) << '\n';
    }
}

void Molrot::print_principal_moments(std::ostream& to) const
{
    chem::Format<char> line;
    line.width(54).fill('-');

    if (!atoms.empty()) {
        to << "\nPrincipal axes and moments of inertia in atomic units:\n"
           << line('-') << '\n';

        chem::Format<double> fix;
        fix.fixed().width(12);

        to << "\t\tA\t     B\t\t  C\n"
           << "Eigenvalue: " << fix(pmom(0)) << ' ' << fix(pmom(1)) << ' '
           << fix(pmom(2)) << '\n';
        to << "     X      ";
        for (std::size_t i = 0; i < pmom.size(); ++i) {
            to << fix(paxis(0, i)) << ' ';
        }
        to << "\n     Y      ";
        for (std::size_t i = 0; i < pmom.size(); ++i) {
            to << fix(paxis(1, i)) << ' ';
        }
        to << "\n     Z      ";
        for (std::size_t i = 0; i < pmom.size(); ++i) {
            to << fix(paxis(2, i)) << ' ';
        }
        to << '\n';
    }
}

void Molrot::print_constants(std::ostream& to)
{
    const double gHz2icm = datum::giga / (datum::c_0 * 100.0);

    chem::Format<char> line;
    line.width(21).fill('-');

    chem::Format<double> fix;
    fix.fixed().width(14);

    arma::vec3 r     = constants();
    std::string symm = symmetry();

    if (r(0) > 0.0) {
        to << "\nRotational constants:\n" << line('-') << '\n';

        if (symm.find("linear") == 0) {
            to << fix(r(0)) << " GHz\n" << fix(r(0) * gHz2icm) << " cm^-1\n\n";
        }
        else {
            double ra = r[0];
            double rb = r[1];
            double rc = r[2];
            to << "\tA\t\tB\t\tC\n"
               << fix(ra) << '\t' << fix(rb) << '\t' << fix(rc) << " GHz\n"
               << fix(ra * gHz2icm) << '\t' << fix(rb * gHz2icm) << '\t'
               << fix(rc * gHz2icm) << " cm^-1\n\n";
        }
        to << "Rotational symmetry number: " << sigma << '\n';
        if (symm.find("atom") == 0) {
            to << "Rotational symmetry: This is an atom\n";
        }
        else {
            to << "Rotational symmetry: " << symm << '\n';
        }
    }
}
