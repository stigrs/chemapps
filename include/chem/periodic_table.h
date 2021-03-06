// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_PERIODIC_TABLE_H
#define CHEM_PERIODIC_TABLE_H

#include <chem/element.h>
#include <stdexcept>
#include <string>

namespace Chem {

// Namespace providing the Periodic Table of Elements.
//
// Source:
//   Zucker, M.A., Kishore, A.R., Sukumar, R., and Dragoset, R.A. (2015),
//   Elemental Data Index (version 2.5). [Online]
//   Available: http://physics.nist.gov/EDI [2016, April 2].
//   National Institute of Standards and Technology, Gaithersburg, MD.
//
namespace Periodic_table {

    struct Bad_atomic_symbol : std::domain_error {
        Bad_atomic_symbol(std::string s) : std::domain_error(s) {}
    };

    Element get_element(const std::string& symbol);
    std::string get_atomic_symbol(const std::string& symbol);
    std::string get_atomic_symbol(int atomic_number);

    int get_max_atomic_number();
    int get_atomic_number(const std::string& symbol);
    int get_mass_number(const std::string& symbol);

    double get_atomic_mass(const std::string& symbol);
    double get_atomic_weight(const std::string& symbol);
    double get_isotope_composition(const std::string& symbol);

    bool atomic_symbol_is_valid(const std::string& symbol);

} // namespace Periodic_table

inline std::string Periodic_table::get_atomic_symbol(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_symbol;
}

inline int Periodic_table::get_atomic_number(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_number;
}

inline int Periodic_table::get_max_atomic_number()
{
    return get_atomic_number("Og");
}

inline int Periodic_table::get_mass_number(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.mass_number;
}

inline double Periodic_table::get_atomic_mass(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_mass;
}

inline double Periodic_table::get_atomic_weight(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_weight;
}

inline double Periodic_table::get_isotope_composition(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.isotope_comp;
}

} // namespace Chem

#endif // CHEM_PERIODIC_TABLE_H
