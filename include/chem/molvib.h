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

#ifndef CHEM_MOLVIB_H
#define CHEM_MOLVIB_H

#include <chem/molrot.h>
#include <srs/array.h>
#include <srs/math.h>
#include <srs/packed.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------

// Error reporting:

struct Molvib_error : std::runtime_error {
    Molvib_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

// Class for handling molecular vibrations.
class Molvib {
public:
    Molvib(Molrot& rot_) : rot(rot_) {}
    Molvib(std::istream& from, const std::string& key, Molrot& rot_) : rot(rot_)
    {
        init(from, key);
    }

    Molvib(const Molvib& vib)
        : hess(vib.hess),
          freqs(vib.freqs),
          mu_freqs(vib.mu_freqs),
          k_fc(vib.k_fc),
          l_cart(vib.l_cart),
          rot(vib.rot)
    {
    }

    // Perform vibrational analysis.
    void analysis(std::ostream& to = std::cout);

    // Get Hessians.
    const srs::packed_dmatrix& get_hessians() const { return hess; }

    // Get vibrational frequencies.
    const srs::dvector& get_freqs() const { return freqs; }

    // Get reduced masses.
    const srs::dvector& get_red_mass() const { return mu_freqs; }

    // Get force constants.
    const srs::dvector& get_force_constant() const { return k_fc; }

    // Get normal coordinates.
    const srs::dcube& get_norm_coords() const { return l_cart; }

    // Get mass-weighted Hessians.
    srs::packed_dmatrix get_mw_hessians() const;

    // Calculate vibrational frequencies from Hessians.
    srs::dvector calc_freqs() const;

    // Calculate Cartesian frequencies from Hessians.
    srs::dvector calc_cart_freqs() const;

    // Calculate zero-point vibrational energy.
    double zero_point_energy() const;

    // Perform normal mode analysis.
    void norm_mode_analysis();

    // Print vibrational frequencies.
    void print(std::ostream& to = std::cout);

private:
    // Initialize.
    void init(std::istream& from, const std::string& key);

    // Set up coordinate vectors for translation and rotation about principal
    // axes of inertia.
    void trans_rot_vec(srs::dcube& dmat, int& n_tr_rot) const;

    // Transform Cartesian Hessians to internal coordinates.
    void trans_hess_int_coord(srs::dcube& dmat,
                              srs::dmatrix& lmat,
                              srs::size_t& n_tr_rot) const;

    // Print Cartesian vibrational frequencies.
    void print_cart_freq(std::ostream& to) const;

    // Print normal modes.
    void print_norm_modes(std::ostream& to) const;

    // Shuffle n_vib orthogogal vectors to the beginning of D matrix.
    void shuffle(srs::dcube& dmat, srs::size_t n_tr_rot) const;

    // Convert frequencies from atomic units to cm^-1.
    void freqs_unit_conv(srs::dvector& vib) const;

    srs::packed_dmatrix hess;  // packed Hessian matrix
    srs::dvector freqs;        // vibrational frequencies
    srs::dvector mu_freqs;     // reduced masses for vibrational modes
    srs::dvector k_fc;         // force constants for vibrational modes
    srs::dcube l_cart;         // Cartesian displacements

    Molrot& rot;
};

#endif  // CHEM_MOLVIB_H
