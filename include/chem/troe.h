////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013-2018 Stig Rune Sellevag. All rights reserved.
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

#ifndef CHEM_TROE_H
#define CHEM_TROE_H

#include <chem/molecule.h>
#include <chem/mol_type.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

// Error reporting:

struct Troe_error : std::runtime_error {
    Troe_error(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Class to perform Troe factorization of strong-collision low-pressure
// limiting rate coefficients for dissociation reactions as presented in
// the papers:
//
//   Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4745-4757.
//   Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758-4775.
//
// The formalism has been extended to include calculation of reduced collision
// integrals using Eq. A4.10 given in
//
//   Forst, W. Unimolecular Reactions; Cambridge University Press, 2003.
//
// Note: Only neutral species can be investigated in this version.
//
class Troe {
public:
    Troe(std::istream& from, Molecule& mol_);

    // Calculate energy dependence factor.
    //
    // Eq. 9.10 in Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758--4775.
    //
    double f_energy(const double temp) const;

    // Calculate anharmonicity factor for m Morse oscillators with
    // dissociation energies equal to the barrier height of reaction.
    //
    // Eq. 5.4 in Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758--4775.
    //
    double f_anharm() const;

    // Calculate rotational factor.
    //
    // Eq. 7.23 (linear molecule, Case II potential), Eq. 7.24 (nonlinear
    // molecule, Case I potential), Eq. 7.26 (linear molecule, Case II
    // potential), and Eq. 7.27 (nonlinear molecule, Case I potential) in
    // Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758--4775.
    //
    double f_rotation(const double temp) const;

    // Calculate free internal rotation factor.
    //
    // Eq. 9.15 in Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758--4775.
    //
    //double f_free_rotor(const double temp) const;

    // Calculate hindered internal rotation factor.
    //
    // Eq. 9.16 in Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758--4775.
    // Hindered internal rotation with E0/V0 <= 3 cannot be treated yet.
    //
    //double f_hind_rotor(const double temp) const;

private:
    Molecule& mol;
    Pot_type pot_type;  // potential type (type 1 or type 2)
    double e_barrier;   // energy barrier for the reaction
    double imom_ratio;  // moment of inertia ratio
    int n_morse_osc;    // number of Morse oscillators

	double zpe; // zero-point vibrational energy
};

inline double Troe::f_anharm() const
{
    double s = mol.get_vib().get_freqs().size();

    return std::pow((s - 1.0) / (s - 1.5), n_morse_osc);
}

#endif  // CHEM_TROE_H
