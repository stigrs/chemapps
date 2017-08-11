/**
   @file mcmm.cpp
 
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

#include <map>
#include <limits>
#include <chem/mcmm.h>
#include <chem/input.h>
#include <chem/utils.h>
#include <chem/mopac.h>


template<class Pot>
Mcmm<Pot>::Mcmm(std::istream& from,
           Molecule& mol_,
           const std::string& key,
           bool verbose_)
    : mol(mol_), verbose(verbose_)
{
    typedef std::map<std::string, Input>::iterator       Input_iter;
    typedef std::map<std::string, Input>::const_iterator Cinput_iter;

    // Read input data:

    const double emin_def = -std::numeric_limits<double>::min();
    
    std::map<std::string, Input> input_data;
    input_data["xtol"]      = Input(xtol, 1.0e-5);
    input_data["etol"]      = Input(etol, 1.0e-8);
    input_data["emin"]      = Input(emin, emin_def);
    input_data["emax"]      = Input(emax, 0.0);
    input_data["rmin"]      = Input(rmin, 0.7414); // experimental r(H-H)
    input_data["temp"]      = Input(temp, 298.15);
    input_data["maxiter"]   = Input(maxiter, 1000);
    input_data["maxreject"] = Input(maxreject, 100);
    input_data["nminima"]   = Input(nminima, 20);
    input_data["seed"]      = Input(seed, 0);
    
    bool found = chem::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            } 
            else {
                Input_iter it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    } 
    else {
        throw Mol_error("cannot find " + key + " section");
    }
    
    // Check if initialized:

    for (Cinput_iter it = input_data.begin(); it != input_data.end(); ++it) {
        if (! it->second.is_init()) {
            throw Mol_error(it->first + " not initialized");
        }
    }

    // Seed the random number engine:

    if (seed == 0) {
        std::random_device rd;
        std::seed_seq seed_seq_{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(seed_seq_);
    }
    else {
        mt.seed(seed);
    }

    std::uniform_real_distribution<> rnd_uni_real(0.0, 1.0);
    for (int i = 0; i < 10; ++i) {
        std::cout << rnd_uni_real(mt) << '\n';
    }
}

template class Mcmm<Mopac>;
