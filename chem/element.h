/**
   @file element.h
 
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

#ifndef CHEM_ELEMENT_H
#define CHEM_ELEMENT_H

#include <string>


/// Struct for holding an element in the Periodic Table of Elements.
struct Element {
    std::string atomic_symbol; ///< atomic or isotopic symbol
    int         atomic_number; ///< atomic number
    int         mass_number;   ///< mass number
    double      atomic_mass;   ///< atomic mass
    double      atomic_weight; ///< atomic weight
    double      isotope_comp;  ///< isotopic composition
};

#endif /* CHEM_ELEMENT_H */

