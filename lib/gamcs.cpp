// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gamcs.h>
#include <chem/gaussian.h>
#include <chem/mopac.h>
#include <chem/io.h>
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <string>
#include <limits>

template <class Pot>
Chem::Gamcs<Pot>::Gamcs(std::istream& from, std::ostream& to) : mol(from, to)
{
    using namespace Stdutils;
    const std::string key = "Gamcs";

    // Read input:

    dist_min = 0.5; // experimental r(H-H) = 0.74 angstrom
    dist_max = 2.2;

    energy_min = -std::numeric_limits<double>::max();
    energy_max = 0.0;
    energy_var = 1.0e-3;
    energy_tol = 1.0e-3;
    ediff_global = 0.0;

    fit_sum_lim = 1.2;
    xyz_rmsd = 0.2;

    prob_cross = 0.95;
    prob_mut = 0.5;

    pop_size = 20;
    max_gen = 200;
    min_gen = 10;
    max_mut_tors = 2;
    cross_trials = 20;
    mut_trials = 100;
    seed = 0;

    select_method = "roulette";

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "dist_min", dist_min, dist_min);
        get_token_value(from, pos, "dist_max", dist_max, dist_max);
        get_token_value(from, pos, "energy_min", energy_min, energy_min);
        get_token_value(from, pos, "energy_max", energy_max, energy_max);
        get_token_value(from, pos, "energy_var", energy_var, energy_var);
        get_token_value(from, pos, "energy_tol", energy_tol, energy_tol);
        get_token_value(from, pos, "fit_sum_lim", fit_sum_lim, fit_sum_lim);
        get_token_value(from, pos, "xyz_rmsd", xyz_rmsd, xyz_rmsd);
        get_token_value(from, pos, "prob_cross", prob_cross, prob_cross);
        get_token_value(from, pos, "prob_mut", prob_mut, prob_mut);
        get_token_value(from, pos, "pop_size", pop_size, pop_size);
        get_token_value(from, pos, "max_gen", max_gen, max_gen);
        get_token_value(from, pos, "min_gen", min_gen, min_gen);
        get_token_value(from, pos, "max_mut_tors", max_mut_tors, max_mut_tors);
        get_token_value(from, pos, "cross_trials", cross_trials, cross_trials);
        get_token_value(from, pos, "mut_trials", mut_trials, mut_trials);
        get_token_value(from, pos, "seed", seed, seed);
        get_token_value(from, pos, "select_method", select_method,
                        select_method);
    }

    // Validate input:

    Assert::dynamic(dist_min >= 0.5, "bad dist_min < 0.5");
    Assert::dynamic(dist_max > dist_min, "bad dist_max <= dist_min");
    Assert::dynamic(energy_var > 0.0, "bad energy_var <= 0.0");
    Assert::dynamic(energy_tol > 0.0, "bad energy_tol <= 0.0");
    Assert::dynamic(fit_sum_lim > 1.0, "bad fit_sum_lim <= 1.0");
    Assert::dynamic(xyz_rmsd > 0.0, "bad xyz_rmsd <= 0.0");
    Assert::dynamic(prob_cross > 0.0 && prob_cross < 1.0,
                    "bad prob_cross != (0, 1)");
    Assert::dynamic(prob_mut > 0.0 && prob_mut < 1.0, "bad prob_mut != (0, 1)");
    Assert::dynamic(pop_size > 1, "bad pop_size <= 1");
    Assert::dynamic(max_gen > pop_size, "bad max_gen <= pop_size");
    Assert::dynamic(min_gen > 1, "bad min_gen <= 1");
    Assert::dynamic(max_mut_tors >= 1, "bad max_mut_tors < 1");
    Assert::dynamic(cross_trials >= 1, "bad cross_trials < 1");
    Assert::dynamic(mut_trials >= 1, "bad mut_trials < 1");
    Assert::dynamic(select_method == "random" || select_method == "roulette" ||
                        select_method == "elite",
                    "bad select_method: " + select_method);

    // Seed the random number engine:

    if (seed == 0) {
        std::random_device rd;
        std::seed_seq seed_seq_{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(seed_seq_);
    }
    else {
        mt.seed(seed); // should only be used for testing purposes
    }

    // Print input parameters:
    print_params(to);

    // Initialize potential:

    pot.init(from);

    // Initialize population:

    init_population(to);
}

template <class Pot>
void Chem::Gamcs<Pot>::solve(std::ostream& to)
{
    Stdutils::Format<char> line;
    line.width(58).fill('-');

    Stdutils::Format<int> ifix;
    ifix.fixed().width(4);

    Stdutils::Format<double> dfix;
    dfix.fixed().width(12).precision(6);

    Stdutils::Format<double> sci;
    sci.scientific().width(8).precision(2);

    to << "Evolution:\n"
       << line('-') << '\n'
       << "Iter  E(best)       E(diff)   E(tol)    Optimization\n"
       << line('-') << std::endl;

    std::uniform_real_distribution<> rnd_uni_real(0.0, 1.0);

    int iter = 0;
    int nsuccess = 0;
    int nfailed = 0;
    bool converged = false;
    std::string status = "";

    while (iter < max_gen && !converged) {
        // Select offsprings:
        Chem::Molecule child1(mol);
        Chem::Molecule child2(mol);
        select_parents(child1, child2);

        // Perform crossover:
        if (rnd_uni_real(mt) < prob_cross) {
            crossover(child1, child2);
        }

        // Perform mutation, making sure that child1 and child2 are sensible
        // and are not in the blacklist:
        if (rnd_uni_real(mt) < prob_mut) {
            mutate(child1);
            mutate(child2);
        }
        blacklist.push_back(
            Chem::Conformer(child1.elec().energy(), child1.get_xyz()));
        blacklist.push_back(
            Chem::Conformer(child2.elec().energy(), child2.get_xyz()));

        // Perform local optimization:
        pot.run(child1);
        pot.run(child2);

        // Update blacklist:
        blacklist.push_back(
            Chem::Conformer(child1.elec().energy(), child1.get_xyz()));
        blacklist.push_back(
            Chem::Conformer(child2.elec().energy(), child2.get_xyz()));

        // Update population if energies are sensible:
        if (child1.elec().energy() >= energy_min &&
            child1.elec().energy() < energy_max &&
            child2.elec().energy() >= energy_min &&
            child2.elec().energy() < energy_max) {
            population.push_back(
                Chem::Conformer(child1.elec().energy(), child1.get_xyz()));
            population.push_back(
                Chem::Conformer(child2.elec().energy(), child2.get_xyz()));

            sort_population();

            population.pop_back(); // delete high-energy candidates
            population.pop_back();

            compute_fitness();
            min_energy.push_back(population[0].energy); // update energy log

            status = "success";
            ++nsuccess;
        }
        else {
            status = "failed";
            ++nfailed;
        }
        // Check convergence:
        if (energy_converged(iter)) {
            converged = true;
        }
        to << ifix(iter + 1) << "  " << dfix(population[0].energy) << "  "
           << sci(ediff_global) << "  " << sci(energy_var) << "  " << status
           << std::endl;

        ++iter;
    }
    to << line('-') << '\n'
       << "Number of successful trials: " << nsuccess << '\n'
       << "Number of failed trials:     " << nfailed << "\n\n";

    line.width(13).fill('-');
    to << "Local minima:\n" << line('-') << '\n';
    print_population(to);
    print_global_minimum(to);
}

//------------------------------------------------------------------------------

template <class Pot>
void Chem::Gamcs<Pot>::init_population(std::ostream& to)
{
    Stdutils::Format<char> line;
    line.width(58).fill('-');

    Stdutils::Format<int> ifix;
    ifix.fixed().width(4);

    Stdutils::Format<double> dfix;
    dfix.fixed().width(12).precision(6);

    to << "Initialization:\n"
       << line('-') << '\n'
       << "Iter  E(curr)       E(best)       Optimization\n"
       << line('-') << std::endl;

    int ipop = 0;
    int iter = 0;
    int nsuccess = 0;
    int nfailed = 0;
    double ecurr = 0.0;
    double ebest = 0.0;

    std::string status = "";

    while (ipop < pop_size && iter < max_gen) {
        // Generate new random structure with sensible geometry:
        Chem::Molecule m(mol);
        gen_rand_conformer(m);
        if (!geom_sensible(m)) {
            continue;
        }

        // Add starting structure for local optimization to blacklist:
        blacklist.push_back(Chem::Conformer(m.elec().energy(), m.get_xyz()));

        // Perform local optimization:
        pot.run(m);
        ecurr = m.elec().energy();

        if (m.elec().energy() >= energy_min && m.elec().energy() < energy_max) {
            // Add optimized structure to blacklist and population:
            blacklist.push_back(
                Chem::Conformer(m.elec().energy(), m.get_xyz()));
            population.push_back(
                Chem::Conformer(m.elec().energy(), m.get_xyz()));
            if (ecurr < ebest) {
                ebest = ecurr;
            }
            status = "success";
            ++nsuccess;
            ++ipop;
        }
        else {
            status = "failed";
            ++nfailed;
        }
        to << ifix(iter + 1) << "  " << dfix(ecurr) << "  " << dfix(ebest)
           << "  " << status << std::endl;
        ++iter;
    }
    to << line('-') << '\n'
       << "Number of successful trials: " << nsuccess << '\n'
       << "Number of failed trials:     " << nfailed << '\n'
       << std::endl;

    sort_population();
    compute_fitness();
    estart = population[0].energy; // save initial energy of global minimum

    line.width(19).fill('-');
    to << "Initial population:\n" << line('-') << '\n';
    print_population(to);
}

template <class Pot>
void Chem::Gamcs<Pot>::gen_rand_conformer(Chem::Molecule& m)
{
    std::uniform_int_distribution<> rnd_uni_int(1, max_mut_tors);
    int n_mut_tors = rnd_uni_int(mt);

    for (int it = 0; it < n_mut_tors; ++it) {
        // Select a random dihedral angle:
        auto moiety = select_rand_dihedral(m);

        // Apply random variation to dihedral angle:
        std::uniform_real_distribution<> rnd_uni_real(-179.0, 180.0);
        double delta = rnd_uni_real(mt);
        m.geom().rotate_moiety(moiety, delta);
    }
}

template <class Pot>
std::vector<Index>
Chem::Gamcs<Pot>::select_rand_dihedral(const Chem::Molecule& m)
{
    auto connect = m.geom().get_connectivities();
    std::uniform_int_distribution<std::size_t> rnd_uni_int(2,
                                                           connect.size() - 1);
    auto index = rnd_uni_int(mt);
    auto dihedral = connect[index];

    std::vector<Index> res;
    for (std::size_t i = 2; i < connect.size(); ++i) {
        if (connect[i] == dihedral) {
            res.push_back(i);
        }
    }
    return res;
}

template <class Pot>
bool Chem::Gamcs<Pot>::geom_sensible(const Chem::Molecule& m) const
{
    bool geom_ok = true;
    Numlib::Mat<double> dist_mat;
    Numlib::pdist_matrix(dist_mat, m.get_xyz());
    for (auto v : dist_mat) {
        if (v > 0.0 && v < dist_min) { // avoid too close atoms
            geom_ok = false;
            break;
        }
    }
    if (geom_ok) { // avoid too long bond distances
        for (std::size_t i = 0; i < m.num_atoms(); ++i) {
            if (m.geom().get_distance(i) >= dist_max) {
                geom_ok = false;
                break;
            }
        }
    }
    return geom_ok;
}

template <class Pot>
bool Chem::Gamcs<Pot>::is_blacklisted(const Numlib::Mat<double>& xyz) const
{
    bool res = false;
    for (const auto& conformer : blacklist) {
        if (Numlib::kabsch_rmsd(conformer.xyz, xyz) < xyz_rmsd) {
            res = true;
            break;
        }
    }
    return res;
}

template <class Pot>
void Chem::Gamcs<Pot>::compute_fitness()
{
    fitness.clear();

    double emin = population[0].energy;
    double emax = population.end()[-1].energy;
    double ediff = std::abs(emax - emin);

    double fi = 0.0;
    for (const auto& p : population) {
        if (ediff < energy_var) {
            fi = 1.0;
        }
        else {
            fi = (emax - p.energy) / ediff;
        }
        fitness.push_back(fi);
    }
}

template <class Pot>
void Chem::Gamcs<Pot>::select_parents_random(std::size_t& indx1,
                                             std::size_t& indx2)
{
    std::uniform_int_distribution<std::size_t> rnd_uni_int(
        0, population.size() - 1);

    indx1 = rnd_uni_int(mt);
    indx2 = 0;

    bool equal = true;
    while (equal) {
        indx2 = rnd_uni_int(mt);
        if (indx1 != indx2) {
            equal = false;
        }
    }
}

template <class Pot>
void Chem::Gamcs<Pot>::select_parents_roulette(std::size_t& indx1,
                                               std::size_t& indx2)
{
    double fit_sum = std::accumulate(fitness.begin(), fitness.end(), 0.0);

    if (fit_sum <= fit_sum_lim) {
        std::uniform_int_distribution<std::size_t> rnd_uni_int(
            1, fitness.size() - 1);
        indx1 = 0;               // select fittest conformer
        indx2 = rnd_uni_int(mt); // and a random conformer
    }
    else {
        indx1 = roulette_select();
        bool equal = true;
        while (equal) {
            indx2 = roulette_select();
            if (indx1 != indx2) {
                equal = false;
            }
        }
    }
}

template <class Pot>
std::size_t Chem::Gamcs<Pot>::roulette_select()
{
    // Algorithm: Fitness proportionate selection (roulette wheel)
    // https://en.wikipedia.org/wiki/Fitness_proportionate_selection

    std::uniform_real_distribution<> rnd_uni_real(0.0, 1.0);
    double rnd = rnd_uni_real(mt);

    double fit_sum = std::accumulate(fitness.begin(), fitness.end(), 0.0);

    std::vector<double> p_ind(fitness.size());
    for (std::size_t i = 0; i < fitness.size(); ++i) {
        p_ind[i] = fitness[i] / fit_sum;
    }
    double cdf = 0.0;
    for (std::size_t i = 0; i < fitness.size(); ++i) {
        cdf += p_ind[i];
        if (rnd < cdf) {
            return i;
        }
    }
    return fitness.size() - 1; // return last item if rounding errors occur
}

template <class Pot>
void Chem::Gamcs<Pot>::crossover(Chem::Molecule& child1, Chem::Molecule& child2)
{
    // Note: Only torsional modes are considered in the crossing-over
    // procedure.

    std::uniform_int_distribution<std::size_t> rnd_uni_int(3,
                                                           mol.num_atoms() - 1);
    std::size_t pt = rnd_uni_int(mt); // crossover point (see Figure S1)

    Chem::Molecule parent1(child1);
    Chem::Molecule parent2(child2);

    bool accepted = false;
    int iter = 0;
    while (!accepted && iter < cross_trials) {
        for (std::size_t i = pt; i < mol.num_atoms(); ++i) {
            child1.geom().set_dihedral(i, parent2.geom().get_dihedral(i));
        }
        if (geom_sensible(child1)) {
            accepted = true;
        }
        ++iter;
    }
    accepted = false;
    iter = 0;
    while (!accepted && iter < cross_trials) {
        for (std::size_t i = pt; i < mol.num_atoms(); ++i) {
            child2.geom().set_dihedral(i, parent1.geom().get_dihedral(i));
        }
        if (geom_sensible(child2)) {
            accepted = true;
        }
        ++iter;
    }
}

template <class Pot>
void Chem::Gamcs<Pot>::mutate(Chem::Molecule& child)
{
    int iter = 0;
    while (iter < mut_trials) {
        gen_rand_conformer(child);
        if (is_blacklisted(child.get_xyz())) {
            continue;
        }
        if (geom_sensible(child)) {
            break;
        }
        ++iter;
    }
}

template <class Pot>
void Chem::Gamcs<Pot>::print_params(std::ostream& to) const
{
    Stdutils::Format<char> line;
    line.width(52).fill('=');

    to << line('=') << '\n'
       << "Genetic Algorithm Molecular Conformer Search (GAMCS)\n"
       << line('=') << "\n\n";

    line.width(52).fill('-');
    to << "Input parameters:\n"
       << line('-') << '\n'
       << "Smallest atom-atom distance allowed:   " << dist_min << '\n'
       << "Largest bond distance allowed:         " << dist_max << '\n'
       << "Geometry RMSD for unique conformers:   " << xyz_rmsd << '\n'
       << "Smallest energy allowed:               " << energy_min << '\n'
       << "Largest energy allowed:                " << energy_max << '\n'
       << "Energy convergence threshold:          " << energy_tol << '\n'
       << "Threshold for variance in energies:    " << energy_var << '\n'
       << "Threshold for sum of fitness values:   " << fit_sum_lim << '\n'
       << "Probability for crossing-over:         " << prob_cross << '\n'
       << "Probability for mutations:             " << prob_mut << '\n'
       << "Population size:                       " << pop_size << '\n'
       << "Maximum number of generations:         " << max_gen << '\n'
       << "Minimum number of generations:         " << min_gen << '\n'
       << "Maximum number of torsional mutations: " << max_mut_tors << '\n'
       << "Maximum number of crossover trials:    " << cross_trials << '\n'
       << "Maximum number of mutation trials:     " << mut_trials << '\n'
       << "Selection method:                      " << select_method << '\n'
       << std::endl;
}

template <class Pot>
void Chem::Gamcs<Pot>::print_population(std::ostream& to) const
{
    Stdutils::Format<double> fix;
    fix.fixed().width(12).precision(6);

    for (std::size_t i = 0; i < population.size(); ++i) {
        to << "Conformer: " << i + 1 << '\n'
           << "Energy: " << fix(population[i].energy) << '\n';
        Chem::print_geometry(to, mol.atoms(), population[i].xyz);
        to << std::endl;
    }
}

template <class Pot>
void Chem::Gamcs<Pot>::print_global_minimum(std::ostream& to) const
{
    Stdutils::Format<char> line;
    line.width(25).fill('-');

    Stdutils::Format<double> fix;
    fix.fixed().width(12).precision(6);

    to << "Estimated global minimum:\n"
       << line('-') << '\n'
       << "E(start): " << fix(estart) << '\n'
       << "E(final): " << fix(population[0].energy) << '\n';
    Chem::print_geometry(to, mol.atoms(), population[0].xyz);
    to << std::endl;
}

template class Chem::Gamcs<Chem::Gaussian>;
template class Chem::Gamcs<Chem::Mopac>;
