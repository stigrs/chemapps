// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_ENERGY_LEVELS_H
#define CHEM_ENERGY_LEVELS_H

#include <vector>

namespace Chem {

namespace Energy_levels {

    // Calculate harmonic oscillator energy levels up to a maximum energy.
    std::vector<double> harmonic_oscillator(double freq, double emax);

    // Calculate free rotor energy levels up to a maximum energy.
    std::vector<double> free_rotor(double rotc, double emax);

    // Calculate hindered rotor energy levels up to a maximum energy.
    //
    // For energies less than 1.5 times the barrier height the hindered rotor
    // energy levels are computed using the method of Barker and Shovlin
    // (Chem. Phys. Lett., 2004, vol. 383, pp. 203-207), while free rotor
    // energy levels are used for higher energies and if the barrier height is
    // less than 1.0 cm^-1.
    //
    // Note: Only implemented for one-dimensional hindered rotors.
    //
    std::vector<double>
    hindered_rotor(double sigma, double rotc, double barrier, double emax);

} // namespace Energy_levels

} // namespace Chem

#endif // CHEM_ENERGY_LEVELS_H

