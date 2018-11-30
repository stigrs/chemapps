// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_ELECTRONIC_H
#define CHEM_ELECTRONIC_H

#include <numlib/matrix.h>
#include <iostream>
#include <string>

namespace Chem {

// Class for handling electronic states.
//
class Electronic {
public:
    Electronic()
        : charge_{0}, spin{1}, energy_{0.0}, so_degen{1}, so_energy{0.0}
    {
    }

    Electronic(std::istream& from, const std::string& key);

    // Copy semantics:
    Electronic(const Electronic&) = default;
    Electronic& operator=(const Electronic&) = default;

    // Move semantics:
    Electronic(Electronic&&) = default;
    Electronic& operator=(Electronic&&) = default;

    ~Electronic() = default;

    // Get properties:

    auto charge() const { return charge_; }
    auto spin_mult() const { return spin; }
    auto energy() const { return energy_; }

    const auto& spin_orbit_degen() const { return so_degen; }
    const auto& spin_orbit_energy() const { return so_energy; }

    // Set properties:

    void set_charge(int value) { charge_ = value; }
    void set_spin_mult(int value) { spin = value; }
    void set_energy(double value) { energy_ = value; }

private:
    int charge_;    // net electronic charge
    int spin;       // spin multiplicity
    double energy_; // electronic energy

    Numlib::Vec<int> so_degen;     // degeneracies of spin-orbit states
    Numlib::Vec<double> so_energy; // energies of spin-orbit states
};

} // namespace Chem

#endif // CHEM_ELECTRONIC_H

