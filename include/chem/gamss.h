// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_GAMSS_H
#define CHEM_GAMSS_H

#include <chem/molecule.h>
#include <chem/conformer.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>

namespace Chem {

// Class providing molecular structure search using a genetic algorithm.
//
template <class Pot>
class Gamss {
public:
    Gamss(std::istream& from, std::ostream& to = std::cout);

    // Run solver.
    void solve(std::ostream& to = std::cout);

private:
    // Initialize population.
    void init_population(std::ostream& to = std::cout);

    // Generate a new random conformer.
    void gen_rand_conformer(Molecule& m);

    // Select random dihedral angle.
    std::vector<int> select_rand_dihedral(const Molecule& m);

    // Check if geometry of random conformer is sensible.
    bool geom_sensible(const Molecule& m) const;

    // Check if geometry is blacklisted.
    bool is_blacklisted(const Numlib::Mat<double>& xyz) const;

    Molecule mol; // molecule
    Pot pot;      // potential function

    double dist_min; // smallest atom-atom distance permitted
    double dist_max; // largest bond distance permitted

    double energy_min; // smallest energy permitted
    double energy_max; // largest energy permitted

    double rmsd_tol_uniq; // geometry tolerance

    int pop_size;     // population size
    int max_iter;     // max number of iterations
    int max_mut_tors; // max number of torsional mutations
    int seed;         // random number generator seed

    std::vector<Conformer> population; // population of optimized structures
    std::vector<Conformer> blacklist;  // population of blacklisted structures

    std::mt19937_64 mt; // random number engine
};

} // namespace Chem

#endif // CHEM_GAMSS_H

