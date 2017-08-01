/**
   @file arma_io.h
   
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

#ifndef CHEM_ARMA_IO_H
#define CHEM_ARMA_IO_H

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>


namespace chem {

    // Error reporting:

    struct Arma_error : std::runtime_error {
        Arma_error(std::string s) : std::runtime_error(s) { }
    };

    // I/O operators for vectors:

    template<typename T>
    void print_vector(std::ostream& to, const arma::Col<T>& a)
    {
        to << a.size() << " [ ";
        for (unsigned i = 0; i < a.size(); ++i) {
            to << a(i) << " ";
            if (! ((i + 1) % 7) && (i != (a.size() - 1))) {
                to << "\n  ";
            }
        }
        to << ']';
    }

    template<typename T>
    void read_vector(std::istream& from, arma::Col<T>& a)
    {
        unsigned n;
        from >> n;
        if (n < 1) {
            throw Arma_error("read_vector: bad size");
        }
        a.resize(n);

        char ch;
        from >> ch;
        if (ch != '[') {
            throw Arma_error("read_vector: '[' missing");
        }
        for (unsigned i = 0; i < n; ++i) {
            from >> a(i);
        }
        from >> ch;
        if (ch != ']') {
            throw Arma_error("read_vector: ']' missing");
        }
    }
} // chem::

#endif /* CHEM_ARMA_IO_H */

