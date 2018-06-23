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

    // Calculate anharmonicity factor for m Morse oscillators with
    // dissociation energies equal to the barrier height of reaction.
    //
    // Eq. 5.4 in Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758--4775.
    //
    double f_anharm() const;

private:
    Molecule& mol;
    Pot_type pot_type;  // potential type
    double e_barrier;
    double imom_ratio;  // moment of inertia ratio
    int n_morse_osc;    // number of Morse oscillators
};

inline double Troe::f_anharm() const
{
    double s = mol.get_vib().get_freqs().size();

    return std::pow((s - 1.0) / (s - 1.5), n_morse_osc);
}

#endif  // CHEM_TROE_H
