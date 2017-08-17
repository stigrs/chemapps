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

#include <chem/input.h>
#include <chem/mcmm.h>
#include <chem/mopac.h>
#include <chem/utils.h>
#include <algorithm>
#include <limits>
#include <map>

template <class Pot>
Mcmm<Pot>::Mcmm(std::istream& from,
                Molecule& mol_,
                const std::string& key,
                bool verbose_)
    : mol(mol_), verbose(verbose_)
{
    // Read input data:

    const double emin_def = -std::numeric_limits<double>::min();

    std::map<std::string, Input> input_data;
    input_data["xtol"]      = Input(xtol, 1.0e-5);
    input_data["etol"]      = Input(etol, 1.0e-8);
    input_data["emin"]      = Input(emin, emin_def);
    input_data["emax"]      = Input(emax, 0.0);
    input_data["rmin"]      = Input(rmin, 0.7414);  // experimental r(H-H)
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
                auto it = input_data.find(token);
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

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
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
        mt.seed(seed);  // should only be used for testing purposes
    }

    std::uniform_real_distribution<> rnd_uni_real(0.0, 1.0);
    for (int i = 0; i < 10; ++i) {
        std::cout << rnd_uni_real(mt) << '\n';
    }
}

template <class Pot>
void Mcmm<Pot>::new_conformer()
{
    do {
        // Select starting geometry:
        arma::mat xnew;
        uniform_usage(xnew);

        // Generate a new conformer:
        // gen_rand_conformer();
    } while (!accept_geom_dist());

#if 0    
    // Perform geometry optimization:
    pot.run(mol);

    // Check acceptance:
    if (accept_energy()) {
        if (! duplicate()) { // store new conformer
            xcurr = mol.get_xyz();
            ecurr = mol.get_elec_energy();
            save_conformer();
            naccept += 1;
        }
    }
#endif
}

template <class Pot>
void Mcmm<Pot>::uniform_usage(arma::mat& xnew)
{
    xnew = xcurr;
    if (!conformers.empty()) {
        std::vector<Conformer> index;
        int istart     = 0;
        int min_nstart = kiter;
        for (unsigned i = 0; i < conformers.size(); ++i) {  // find least used
            if (conformers[i].iter <= min_nstart) {
                index.push_back(conformers[i]);
            }
        }
        if (!conformers.empty()) {  // find least used with lowest energy
            auto res = std::max_element(conformers.begin(), conformers.end());
            int it   = std::distance(conformers.begin(), res);
            double emin_ = conformers[it].energy;
            for (unsigned i = 0; i < index.size(); ++i) {
                if (index[i].energy < emin_) {
                    emin_  = index[i].energy;
                    istart = i;
                }
            }
        }
        conformers[istart].iter += 1;
        xnew = conformers[istart].xyz;
    }
}

template <class Pot>
void Mcmm<Pot>::gen_rand_conformer(Molecule& m) const
{
    // Select a random dihedral angle:

}
template class Mcmm<Mopac>;
