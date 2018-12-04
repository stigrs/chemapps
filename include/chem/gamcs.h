// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_GAMCS_H
#define CHEM_GAMCS_H

#include <chem/molecule.h>
#include <chem/conformer.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>

namespace Chem {

// Class providing molecular structure search using a genetic algorithm.
//
template <class Pot>
class Gamcs {
public:
    Gamcs(std::istream& from, std::ostream& to = std::cout);

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

    // Check if energy of global minimum has converged.
    bool energy_converged(int iter);

    // Compute fitness for a sorted population.
    void compute_fitness();

    // Select parents from population.
    void select_parents(Molecule& parent1, Molecule& parent2);

    // Select parents with uniform probability.
    void select_parents_random(std::size_t& indx1, std::size_t& indx2);

    // Select parents using roulette wheel method.
    void select_parents_roulette(std::size_t& indx1, std::size_t& indx2);

    // Select parents using elite method.
    void select_parents_elite(std::size_t& indx1, std::size_t& indx2);

    // Select index using roulette wheel method.
    std::size_t roulette_select();

    // Perform crossover.
    void crossover(Molecule& child1, Molecule& child2);

    // Perform mutation.
    void mutate(Molecule& child);

    // Sort population.
    void sort_population();

    // Print input parameters.
    void print_params(std::ostream& to) const;

    // Print population.
    void print_population(std::ostream& to) const;

    // Print estimated global minimum.
    void print_global_minimum(std::ostream& to) const;

    Molecule mol; // molecule
    Pot pot;      // potential function

    double dist_min; // smallest atom-atom distance permitted
    double dist_max; // largest bond distance permitted
    double xyz_rmsd; // geometry tolerance

    double energy_min;   // smallest energy permitted
    double energy_max;   // largest energy permitted
    double energy_var;   // lowest energy variance permitted
    double energy_tol;   // energy convergence tolerance
    double ediff_global; // energy difference for global minimum

    double fit_sum_lim; // threshold for sum of fitness values
    double prob_cross;  // probability for crossing
    double prob_mut;    // probability for mutation

    int pop_size;     // population size
    int max_gen;      // max number of generations
    int min_gen;      // min number of generations
    int max_mut_tors; // max number of torsional mutations
    int cross_trials; // max number of crossing trials
    int mut_trials;   // max number of mutation trials
    int seed;         // random number generator seed

    std::string select_method; // selection algorithm

    std::vector<Conformer> population; // population of optimized structures
    std::vector<Conformer> blacklist;  // population of blacklisted structures
    std::vector<double> fitness;       // fitness of population
    std::vector<double> min_energy;    // energies of most stable conformer

    std::mt19937_64 mt; // random number engine
};

template <class Pot>
inline bool Chem::Gamcs<Pot>::energy_converged(int iter)
{
    bool converged = false;
    if (iter > min_gen) {
        std::sort(min_energy.begin(), min_energy.end());
        if (min_energy.size() > min_gen) {
            min_energy.pop_back();
        }
        if (min_energy[0] < energy_min) {
            converged = true;
        }
        double ei = min_energy.end()[-1];
        double e0 = min_energy[0];
        ediff_global = std::abs(ei - e0);
        if (ediff_global < energy_tol) {
            converged = true;
        }
    }
    return converged;
}

template <class Pot>
inline void Gamcs<Pot>::select_parents(Molecule& parent1, Molecule& parent2)
{
    std::size_t indx1 = 0;
    std::size_t indx2 = 0;

    if (select_method == "roulette") {
        select_parents_roulette(indx1, indx2);
    }
    else if (select_method == "elite") {
        select_parents_elite(indx1, indx2);
    }
    else {
        select_parents_random(indx1, indx2);
    }
    parent1.set_xyz(population[indx1].xyz);
    parent2.set_xyz(population[indx2].xyz);
}

template <class Pot>
inline void Gamcs<Pot>::select_parents_elite(std::size_t& indx1,
                                             std::size_t& indx2)
{
    indx1 = 0; // always select the most fit conformers
    indx2 = 1; // (generally not recommended)
}

template <class Pot>
inline void Gamcs<Pot>::sort_population()
{
    std::sort(population.begin(), population.end());
}

} // namespace Chem

#endif // CHEM_GAMCS_H

