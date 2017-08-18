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
        throw Mcmm_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Mcmm_error(it->first + " not initialized");
        }
    }

    // Initialize iterators:

    kiter   = 0;
    nreject = 0;
    naccept = 0;

    // Initialize storage containers:

    xcurr   = mol.get_xyz();
    xglobal = xcurr;
    ecurr   = mol.get_elec_energy();
    eglobal = ecurr;

    // Seed the random number engine:

    if (seed == 0) {
        std::random_device rd;
        std::seed_seq seed_seq_{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(seed_seq_);
    }
    else {
        mt.seed(seed);  // should only be used for testing purposes
    }
}

template <class Pot>
void Mcmm<Pot>::solve()
{
    bool finished = false;
    while (!finished) {
        new_conformer();
        update();
        if (check_exit()) {
            finished = true;
        }
    }
}

template <class Pot>
void Mcmm<Pot>::new_conformer()
{
    Molecule m(mol);
    // Generate a new random conformer by using the uniform usage scheme:
    do {
        arma::mat xnew = m.get_xyz();
        uniform_usage(xnew);
        m.set_xyz(xnew);
        gen_rand_conformer(m);
    } while (!accept_geom_dist(m));  // check geometry constraints

    // Perform geometry optimization:
    pot.run(m);
#if 0    
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
void Mcmm<Pot>::update()
{
    kiter += 1;
}

template <class Pot>
bool Mcmm<Pot>::check_exit() const
{
    bool finished = false;
    if (kiter >= maxiter) {
        finished = true;
    }
    return finished;
}


template <class Pot>
void Mcmm<Pot>::uniform_usage(arma::mat& xnew)
{
    xnew = xcurr;
#if 0
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
#endif
}

template <class Pot>
std::vector<int> Mcmm<Pot>::select_rand_dihedral(const Molecule& m)
{
    std::vector<arma::ivec> connect = m.get_zmat()->get_connectivities();
    std::uniform_int_distribution<> rnd_uni_int(2, connect.size() - 1);
    int index           = rnd_uni_int(mt);
    arma::ivec dihedral = connect[index];

    std::vector<int> res(0);
    for (std::size_t i = 2; i < connect.size(); ++i) {
        if (arma::all((connect[i] == dihedral) == 1)) {
            res.push_back(i);
        }
    }
    return res;
}

template class Mcmm<Mopac>;
