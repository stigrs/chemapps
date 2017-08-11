/**
   @file conformer.h
 
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

#ifndef CHEM_CONFORMER_H
#define CHEM_CONFORMER_H

#include <vector>
#include <armadillo>
#include <chem/element.h>


/// Simple structure for storing conformers.
struct Conformer {
    Conformer(double e, const arma::mat& x) : energy(e), atoms(0), xyz(x) { }

    Conformer(double e, const std::vector<Element>& at, const arma::mat& x)
        : energy(e), atoms(at), xyz(x)
    { }
    
    /// Sort conformers by energy in ascending order.
    bool operator<(const Conformer& c) const { return energy < c.energy; }
    
    /// Sort conformers by energy in descending order.
    bool operator>(const Conformer& c) const { return energy > c.energy; }
    
    double energy;              ///< conformer energy
    std::vector<Element> atoms; ///< atoms
    arma::mat xyz;              ///< Cartesian coordinates
};
    
#endif /* CHEM_CONFORMER_H */
