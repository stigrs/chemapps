// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_UNITS_H
#define CHEM_UNITS_H

#include <iostream>
#include <string>

namespace Chem {

// Namespace providing methods for handling units.
namespace Units {

    enum Type {
        kJ_mol,   // kJ/mol
        kcal_mol, // kcal/mol
        icm,      // cm**-1
        kelvin,
        hartree,
        hertz,
        eV,
        amu,
        kg,
        au
    };

    Type lexer(const std::string& unit);
    void print(std::ostream& to = std::cout);

} // namespace Units

} // namespace Chem

#endif // CHEM_UNITS_H

