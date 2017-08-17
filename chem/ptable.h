/**
   @file ptable.h

   This file is part of ChemApps - A C++ Chemistry Toolkit

   Copyright (C) 2016-2017  Stig Rune Sellevag

   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CHEM_PTABLE_H
#define CHEM_PTABLE_H

#include <stdexcept>
#include <string>

#include <chem/element.h>

/**
   Namespace providing the Periodic Table of Elements.

   Source: http://www.nist.gov/plm/data/comp.cfm
   Downloaded: 2 April 2016
*/
namespace ptable {

struct Bad_atomic_symbol : std::domain_error {
    Bad_atomic_symbol(std::string s) : std::domain_error(s) {}
};

Element get_element(const std::string& symbol);
std::string get_atomic_symbol(const std::string& symbol);

int get_atomic_number(const std::string& symbol);
int get_mass_number(const std::string& symbol);

double get_atomic_mass(const std::string& symbol);
double get_atomic_weight(const std::string& symbol);
double get_isotope_composition(const std::string& symbol);

bool atomic_symbol_is_valid(const std::string& symbol);

}  // ptable::

inline std::string ptable::get_atomic_symbol(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_symbol;
}

inline int ptable::get_atomic_number(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_number;
}

inline int ptable::get_mass_number(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.mass_number;
}

inline double ptable::get_atomic_mass(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_mass;
}

inline double ptable::get_atomic_weight(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_weight;
}

inline double ptable::get_isotope_composition(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.isotope_comp;
}

#endif /* CHEM_PTABLE_H */
