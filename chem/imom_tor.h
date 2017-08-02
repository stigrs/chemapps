/**
   @file imom_tor.h
   
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

#ifndef CHEM_IMOM_TOR_H
#define CHEM_IMOM_TOR_H

#include <iostream>
#include <string>
#include <armadillo>


/**
   Class for calculating reduced moment of inertia for torional modes in
   molecules using the C scheme.

   Algorithm:
   ----------
   The reduced moment of inertia of a symmetrical or unsymmetrical rotating
   top attached to a rigid frame is calculated according to eq 1 in
   Pitzer, K. S. J. Chem. Phys. 1946, vol. 14, pp. 239-243, also known as
   the curvilinear (C) scheme.

   Note:
   -----
   Atoms specifying the rotational axis must not be included in the list
   of atoms specifying the rotating top.
*/
class Imom_tor {
public:
    Imom_tor() { }

    Imom_tor(std::istream& from, const std::string& key);

    ~Imom_tor() { }

private:
    arma::ivec rot_axis = arma::zeros<arma::ivec>(2);
    arma::ivec rot_top = arma::zeros<arma::ivec>(0);
};

#endif /* CHEM_IMOM_TOR_H */
