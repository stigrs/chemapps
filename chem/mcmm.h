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

#include <stdexcept>
#include <string>
#include <iostream>
#include <vector>
#include <armadillo>
#include <chem/molecule.h>

//-----------------------------------------------------------------------------

// Error reporting:

struct Mcmm_error : std::runtime_error {
    Mcmm_error(std::string s) : std::runtime_error(s) { }
};

//-----------------------------------------------------------------------------

/// Class providing Monte Carlo Multiple Minima (MCMM) solver.
template<class Pot>
class Mcmm {
public:
    Mcmm(std::istream& from, 
         Molecule& mol_, 
         const std::string& key = "Mcmm", 
         bool verbose = false);

    ~Mcmm() { }

    /// MCMM solver.
    void solve();

private:
#if 0
    /// Check acceptance of energy.
    bool accept_energy(double enew);

    /// Check if MCMM solver is finished.
    bool check_exit() const;

    /// Update MCMM solver.
    void update();

    /**
       Function for selecting starting geometry by using the uniform 
       usage scheme.
    */
    void uniform_usage(arma::mat& xnew);

    /// Generate a new molecular conformer.
    void new_conformer();

    /// Sort local minima in descending order and cut to nminima.
    void sort_minima();

    /// Save local minima.
    void save_conformer(double energy, const Molecule& mol);
#endif
    Molecule& mol; ///< molecule

    double xtol; ///< absolute error in geometry
    double etol; ///< absolute error in energy
    double emin; ///< lowest energy permitted
    double emax; ///< highest energy permitted
    double rmin; ///< smallest bond distance permitted
    double temp; ///< temperature

    int kiter;     ///< iteration parameter
    int maxiter;   ///< maximum number of iterations (k)
    int maxreject; ///< maximum number of consecutive rejected trials
    int nreject;   ///< iterator for number of consecutive rejected trials
    int naccept;   ///< iterator for number of accepted trials
    int nminima;   ///< number of local minima stored
    int seed;      ///< seed for random number generator

    arma::mat xcurr;   ///< current geometry
    arma::mat xglobal; ///< geometry of global minimum

    double ecurr;   ///< current energy
    double eglobal; ///< energy of global minimum

    std::vector<double>    eminima; ///< array with local energy minima
    std::vector<arma::mat> xminima; ///< array with local geometry minima
}; 

#endif /* CHEM_MCMM_H */

