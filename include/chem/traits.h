// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_TRAITS_H
#define CHEM_TRAITS_H

#include <string>

namespace Chem {

// Enumeration of molecular structure types.
enum Mol_type { atom, linear, nonlinear };

// Enumeration of potential types.
enum Pot_type { type1, type2 };

// Struct for holding molecular formula.
struct Mol_formula {
    std::string atom; // atomic or isotopic symbol
    int stoich;       // stoichiometry
};

}  // namespace Chem

#endif  // CHEM_TRAITS_H

