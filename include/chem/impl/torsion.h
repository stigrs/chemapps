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

#ifndef CHEM_TORSION_H
#define CHEM_TORSION_H

#include <chem/impl/geometry.h>
#include <chem/impl/rotation.h>
#include <numlib/traits.h>
#include <numlib/matrix.h>
#include <iostream>
#include <string>

namespace Chem {

namespace Impl {

    // Class for handling torsional modes in molecules using the CT-Cw scheme.
    //
    // Algorithm:
    // ----------
    // The reduced moment of inertia of a symmetrical or unsymmetrical rotating
    // top attached to a rigid frame is calculated according to eq 1 in
    // Pitzer, K. S. J. Chem. Phys. 1946, vol. 14, pp. 239-243, also known as
    // the curvilinear (C) scheme.
    //
    // The CT-Cw scheme is described in the following paper: Chuang, Y.-Y.;
    //  Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.
    //
    // Note:
    // -----
    // Atoms specifying the rotational axis must not be included in the list
    // of atoms specifying the rotating top.
    //
    class Torsion {
    public:
        Torsion() = delete;

        Torsion(Geometry& g, Rotation& r)
            : geom(g), rot(r), perform_analysis{false}
        {
            alpha = Numlib::zeros<Numlib::Mat<double>>(3, 3);
            top_origo = Numlib::zeros<Numlib::Vec<double>>(3);
            top_com = Numlib::zeros<Numlib::Vec<double>>(3);
        }

        Torsion(std::istream& from,
                const std::string& key,
                Geometry& g,
                Rotation& r);

        // Copy semantics:
        Torsion(const Torsion&) = default;
        Torsion& operator=(const Torsion&) = default;

        // Move semantics:
        Torsion(Torsion&&) = default;
        Torsion& operator=(Torsion&&) = default;

        ~Torsion() = default;

        // Perform torsional mode analysis.
        void analysis(std::ostream& to) const;

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

        Geometry& geom;
        Rotation& rot;

        bool perform_analysis;

        Numlib::Mat<double> alpha; // direction cosines
        Numlib::Mat<double> xyz;   // local copy of Cartesian coordinates

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
    };

    inline double Torsion::symmetry_number() const
    {
        // Eq. 8 in Chuang and Truhlar (2000):
        return tot_minima() / narrow_cast<double>(sigma_tor.size());
    }

} // namespace Impl

} // namespace Chem

#endif // CHEM_TORSION_H
