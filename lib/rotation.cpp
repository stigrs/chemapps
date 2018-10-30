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

#include <chem/impl/rotation.h>
#include <chem/impl/io_support.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <cassert>
#include <cmath>

Chem::Impl::Rotation::Rotation(std::istream& from,
                               const std::string& key,
                               Geometry& g)
    : geom(g), pmom(3), paxis(3, 3)
{
    using namespace Stdutils;

    aligned = false;

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "sigma", sigma, 1);
    }
}
#if 0
void Chem::Impl::Rotation::analysis(std::ostream& to)
{
    if (geom.atoms().size() > 1) {
        rotate_to_principal_axes();
        to << "\nGeometry in principal axes coordinate system:\n";
        Chem::Impl::print_geometry(to, geom.atoms(), geom.cart_coord());
        Chem::Impl::print_center_of_mass(to, center_of_mass());
        Chem::Impl::print_principal_moments(to, pmom, paxis);
        Chem::Impl::print_rot_constants(to, sigma, symmetry(), constants());
    }
}
#endif
Numlib::Vec<double> Chem::Impl::Rotation::constants()
{
    using namespace Numlib::Constants;

    if (!aligned) {
        rotate_to_principal_axes();
    }
    // const double tol = 1.0e-3;
    // const double factor = h_bar / (4.0 * pi * giga * m_u * a_0 * a_0
    // * 1.0e-20);

    Numlib::Vec<double> res(3);
#if 0
    if (geom.atoms().size() > 1) {
        if (std::abs(pmom(0)) < tol) {
            res(0) = factor / pmom(2);
        }
        else {
            res(0) = factor / pmom(0);
            res(1) = factor / pmom(1);
            res(2) = factor / pmom(2);
        }
    }
#endif
    return res;
}
#if 0
std::string Chem::Impl::Rotation::symmetry()
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
        if (geom.atoms().size() == 1) {
            symm = "atom";
        }
    }
    else if (bc) {
        symm = "prolate symmetric top";
        if (geom.atoms().size() == 2) {
            symm = "linear " + symm;
        }
    }
    else if (ab) {
        symm = "oblate symmetric top";
    }
    return symm;
}
#endif
Numlib::Vec<double> Chem::Impl::Rotation::center_of_mass() const
{
    Numlib::Vec<double> com(3);

    for (Index j = 0; j < geom.cart_coord().cols(); ++j) {
        double sum = 0.0;
        for (Index i = 0; i < geom.cart_coord().rows(); ++i) {
            sum += geom.atoms()[i].atomic_mass * geom.cart_coord()(i, j);
        }
        com(j) = sum / geom.tot_mass();
    }
    return com;
}

void Chem::Impl::Rotation::calc_principal_moments()
{
    // Work on a local copy of the Cartesian coordinates:
    Numlib::Mat<double> xyz = geom.cart_coord();

    // Convert geometry units to bohr:
    xyz /= Numlib::Constants::a_0;

    // Move geometry to center of mass:
    auto com = center_of_mass();
    Numlib::translate(xyz, -com(0), -com(1), -com(2));
    std::cout << xyz << std::endl;

    // Compute principal moments:

    paxis = 0.0;
    pmom = 0.0;
    std::cout << paxis << std::endl;
    std::cout << pmom << std::endl;

    if (geom.atoms().size() > 1) {
        for (std::size_t i = 0; i < geom.atoms().size(); ++i) {
            double m = geom.atoms()[i].atomic_mass;
            double x = xyz(i, 0);
            double y = xyz(i, 1);
            double z = xyz(i, 2);
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

