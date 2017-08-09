/**
   @file molvib.cpp
   
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

#include <chem/molvib.h>
#include <chem/utils.h>
#include <chem/datum.h>
#include <chem/arma_io.h>


Molvib::Molvib(std::istream& from, const std::string& key)
{
    bool found = chem::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "vibr") {
                chem::read_vector(from, freqs);
            }
        }
    }
    else {
        throw Molvib_error("cannot find " + key + " section");
    }
}

void Molvib::print(std::ostream& to)
{
    chem::Format<char> line;
    line.width(26).fill('-');

    chem::Format<double> fix;
    fix.fixed().width(8).precision(2);

    if (freqs.size() > 0) {
        int it = 0;
        to << "\nVibrational modes (cm^-1):\n" << line('-') << '\n';
        for (arma::uword i = 0; i < freqs.size(); ++i) {
            to << fix(freqs(i));
            if ((it == 8) && (freqs.size() > 9)) {
                to << '\n';
            }
            it += 1;
        }
        double zpe = zero_point_energy();
        to << "\n\nZero-point vibrational energy: "
           << zpe / datum::au_to_icm << " Hartree\n";
    }
}

