/**
   @file mcmm.h

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

#ifndef CHEM_MCMM_H
#define CHEM_MCMM_H

#include <chem/conformer.h>
#include <chem/math.h>
#include <chem/molecule.h>
#include <armadillo>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------

// Error reporting:

struct Mcmm_error : std::runtime_error {
    Mcmm_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

/// Class providing Monte Carlo Multiple Minima (MCMM) solver.
template <class Pot>
class Mcmm {
public:
    Mcmm(std::istream& from,
         Molecule& mol_,
         const std::string& key = "Mcmm",
         bool verbose_          = false);

    ~Mcmm() {}

    /// MCMM solver.
    void solve();

private:
    /// Check if MCMM solver is finished.
    bool check_exit() const;

    /// Generate a new molecular conformer.
    void new_conformer();

    /// Check if random conformer is ok.
    bool accept_geom_dist() const;

    /**
       Function for selecting starting geometry by using the uniform
       usage scheme.
    */
    void uniform_usage(arma::mat& xnew);

    /// Generate a new random conformer.
    void gen_rand_conformer(Molecule& m) const;

#if 0
    /// Check acceptance of energy.
    bool accept_energy(double enew);

    /// Update MCMM solver.
    void update();

    /// Sort local minima in descending order and cut to nminima.
    void sort_minima();

    /// Save local minima.
    void save_conformer(double energy, const Molecule& mol);
#endif
    Molecule& mol;  ///< molecule
    Pot pot;        ///< potential function

    double xtol;  ///< absolute error in geometry
    double etol;  ///< absolute error in energy
    double emin;  ///< lowest energy permitted
    double emax;  ///< highest energy permitted
    double rmin;  ///< smallest bond distance permitted
    double temp;  ///< temperature

    int maxiter;    ///< maximum number of iterations (k)
    int maxreject;  ///< maximum number of consecutive rejected trials
    int nminima;    ///< number of local minima stored
    int seed;       ///< random number generator seed

    int kiter;    ///< iteration parameter
    int nreject;  ///< iterator for number of consecutive rejected trials
    int naccept;  ///< iterator for number of accepted trials

    arma::mat xcurr;    ///< current geometry
    arma::mat xglobal;  ///< geometry of global minimum

    double ecurr;    ///< current energy
    double eglobal;  ///< energy of global minimum

    std::vector<Conformer> conformers;  ///< array with local energy minima

    bool verbose;

    std::mt19937_64 mt;  ///< random number engine
};

template <class Pot>
inline bool Mcmm<Pot>::accept_geom_dist() const
{
    bool geom_ok = true;
    arma::mat dist_mat;
    chem::pdist_matrix(dist_mat, mol.get_xyz());
    if (dist_mat.min() < rmin) {  // avoid too close atoms
        geom_ok = false;
    }
    return geom_ok;
}

#endif /* CHEM_MCMM_H */
