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

#include <chem/molrot.h>
#include <srs/array.h>
#include <srs/math.h>
#include <iostream>
#include <stdexcept>
#include <string>

// Error reporting:

struct Torsion_error : std::runtime_error {
    Torsion_error(std::string s) : std::runtime_error(s) {}
};

//
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
    Torsion(Molrot& rot_) : rot(rot_) {}
    Torsion(std::istream& from, const std::string& key, Molrot& rot_);

    Torsion(const Torsion& tor);

    // Perform torsional mode analysis.
    void analysis(std::ostream& to = std::cout);

    // Get total number of minima (eq 1 in C&T, 2000).
    int tot_minima() const { return srs::sum(sigma_tor); }

    // Calculate effective symmetry number.
    double symmetry_number() const;

    // Calculate reduced moment of inertia.
    double red_moment_of_inertia();

    // Calculate effective moment of inertia.
    double eff_moment_of_inertia() const;

    // Calculate rotational constant for torsional mode.
    srs::dvector constant() const;

    srs::dvector get_pot_coeff() const { return pot_tor; }
    srs::dvector get_freqs() const { return freq_tor; }

private:
    // Initialize input data.
    void init(std::istream& from, const std::string& key);

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

    srs::dmatrix xyz_;              // local copy of xyz
    srs::dmatrix alpha(3, 3, 0.0);  // direction cosines

    srs::ivector rot_axis(2, 0);  // rotational axis
    srs::ivector rot_top;         // rotating top moiety
    srs::ivector sigma_tor;       // symmetry number

    srs::dvector rmi_tor;   // red. moment of inertia
    srs::dvector pot_tor;   // potential coefficients
    srs::dvector freq_tor;  // torsional frequencies

    srs::dvector x_axis(3, 0.0);  // x axis of rot. top
    srs::dvector y_axis(3, 0.0);  // y axis of rot. top
    srs::dvector z_axis(3, 0.0);  // z axis of rot. top

    srs::dvector top_origo(3, 0.0);  // origo
    srs::dvector top_com(3, 0.0);    // center of mass

    double am;  // moment of inertia of rotating top
    double bm;  // xz product of inertia
    double cm;  // yz product of inertia
    double um;  // off-balance factor

    bool perform_torsional_analysis = false;

    Molrot& rot;
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
    return tot_minima() / static_cast<double>(sigma_tor.size());
}

#endif  // CHEM_TORSION_H
