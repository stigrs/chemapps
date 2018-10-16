// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_ELEMENT_H
#define CHEM_ELEMENT_H

#include <string>

namespace Chem {

// Struct for holding an element in the Periodic Table of Elements.
struct Element {
    std::string atomic_symbol; // atomic or isotopic symbol
    int atomic_number;         // atomic number
    int mass_number;           // mass number
    double atomic_mass;        // atomic mass
    double atomic_weight;      // atomic weight
    double isotope_comp;       // isotopic composition
};

}  // namespace Chem

#endif  // CHEM_ELEMENT_H

