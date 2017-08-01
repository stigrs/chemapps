/**
   @file molvib.h
   
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

#ifndef CHEM_MOLVIB_H
#define CHEM_MOLVIB_H

#include <iostream>
#include <stdexcept>
#include <string>
#include <armadillo>

//-----------------------------------------------------------------------------

// Error reporting:

struct Molvib_error : std::runtime_error {
    Molvib_error(std::string s) : std::runtime_error(s) { }
};

//-----------------------------------------------------------------------------

/// Class for handling molecular vibrations.
class Molvib {
public:
    Molvib() { }

    Molvib(std::istream& from, const std::string& key);

    ~Molvib() { }

    /// Get vibrational frequencies.
    const arma::vec& get_freqs() const { return freqs; }

    /// Calculate zero-point vibrational energy.
    double zero_point_energy() const { return 0.5 * arma::sum(freqs); }

    /// Print vibrational frequencies.
    void print(std::ostream& to = std::cout);

private:
    arma::vec freqs = arma::zeros<arma::vec>(0);
};

#endif /* CHEM_MOLVIB_H */
