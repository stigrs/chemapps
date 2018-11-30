// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gaussian.h>
#include <chem/mcmm.h>
#include <chem/io.h>
#include <chem/mopac.h>
#include <stdutils/stdutils.h>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>

template <class Pot>
Chem::Mcmm<Pot>::Mcmm(std::istream& from,
                      const Chem::Molecule& mol_,
                      const std::string& key,
                      bool verbose_)
    : mol(mol_), verbose(verbose_)
{
    // Read input data:

    using namespace Stdutils;

    const double emin_def = -std::numeric_limits<double>::max();

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "xtol", xtol, 5.0e-2);
        get_token_value(from, pos, "etol", etol, 1.0e-2);
        get_token_value(from, pos, "emin", emin, emin_def);
        get_token_value(from, pos, "emax", emax, 0.0);
        get_token_value(from, pos, "rmin", rmin, 0.5);
        get_token_value(from, pos, "temp", temp, 298.15);
        get_token_value(from, pos, "maxiter", maxiter, 500u);
        get_token_value(from, pos, "miniter", miniter, 50u);
        get_token_value(from, pos, "maxreject", maxreject, 100u);
        get_token_value(from, pos, "nminima", nminima, 20u);
        get_token_value(from, pos, "seed", seed, 0);
    }

    // Validate input:

    Assert::dynamic(xtol > 0.0, "bad xtol <= 0.0");
    Assert::dynamic(etol > 0.0, "bad etol <= 0.0");
    Assert::dynamic(rmin > 0.0, "bad rmin <= 0.0");
    Assert::dynamic(temp > 0.0, "bad temp <= 0.0");
    Assert::dynamic(maxiter >= nminima, "bad maxiter < nminima");
    Assert::dynamic(maxiter >= miniter, "bad maxiter < miniter");
    Assert::dynamic(maxreject >= 1, "bad maxreject < 1");
    Assert::dynamic(nminima >= 1, "bad nminima < 1");

    // Initialize potential:

    pot.init(from);

    // Initialize iterators:

    kiter = 0;
    nreject = 0;
    naccept = 0;

    // Initialize storage containers:

    xcurr = mol.get_xyz();
    ecurr = mol.elec().energy();

    // Seed the random number engine:

    if (seed == 0) {
        std::random_device rd;
        std::seed_seq seed_seq_{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(seed_seq_);
    }
    else {
        mt.seed(seed); // should only be used for testing purposes
    }
}

template <class Pot>
void Chem::Mcmm<Pot>::solve(std::ostream& to)
{
    global_min_found = false;
    while (!global_min_found) {
        new_conformer();
        update();
        if (check_exit()) {
            sort_conformers();
            global_min_found = true;
        }
    }
    if (verbose) {
        double eglobal_min = *std::min_element(eglobal.begin(), eglobal.end());
        Stdutils::Format<char> line;
        line.width(41).fill('=');
        Stdutils::Format<int> ifix;
        ifix.fixed().width(5);
        Stdutils::Format<double> dfix;
        dfix.fixed().width(8).precision(2);
        to << "Monte Carlo Multiple Minima (MCMM) Solver\n"
           << line('=') << '\n';
        to << "Temperature:\t" << dfix(temp) << '\n'
           << "Iterations:\t" << ifix(kiter) << " out of " << maxiter << '\n'
           << "Rejections:\t" << ifix(nreject) << " out of " << maxreject
           << "\n\n";
        dfix.fixed().width(12).precision(6);
        line.width(15).fill('-');
        to << "Global minimum:\n"
           << line('-') << '\n'
           << "Energy: " << dfix(eglobal_min) << '\n';
        Chem::print_geometry(to, mol.atoms(), xglobal);
        to << '\n';

        line.width(13).fill('-');
        to << "Local minima:\n" << line('-') << '\n';
        for (std::size_t i = 0; i < conformers.size(); ++i) {
            to << "Conformer: " << i + 1 << '\n'
               << "Energy: " << dfix(conformers[i].energy) << '\n';
            Chem::print_geometry(to, mol.atoms(), conformers[i].xyz);
            to << '\n';
        }
    }
}

template <class Pot>
bool Chem::Mcmm<Pot>::check_exit() const
{
    bool finished = false;
    if (ecurr < emin) {
        finished = true;
    }
    if (kiter >= maxiter) {
        finished = true;
    }
    if (nreject >= maxreject) {
        finished = true;
    }
    if (eglobal.size() > 1) {
        double ediff = eglobal.end()[-1] - eglobal.end()[-2];
        if ((std::abs(ediff) < etol) && (kiter >= miniter)) {
            finished = true;
        }
    }
    return finished;
}

template <class Pot>
bool Chem::Mcmm<Pot>::accept_energy(double enew)
{
    bool accept = false;
    double ediff = enew - ecurr;
    if (enew > emax) {
        nreject += 1;
    }
    else if (ediff < 0.0) {
        accept = true;
    }
    else {
        double h = std::exp(-ediff / temp);
        std::uniform_real_distribution<> rnd_real_uni(0.0, 1.0);
        if (h > rnd_real_uni(mt)) {
            accept = true;
        }
        else {
            nreject += 1;
        }
    }
    return accept;
}

template <class Pot>
bool Chem::Mcmm<Pot>::accept_geom_dist(const Chem::Molecule& m) const
{
    bool geom_ok = true;
    Numlib::Mat<double> dist_mat;
    Numlib::pdist_matrix(dist_mat, m.get_xyz());
    for (auto v : dist_mat) {
        if (v > 0.0 && v < rmin) { // avoid too close atoms
            geom_ok = false;
        }
    }
    return geom_ok;
}

template <class Pot>
bool Chem::Mcmm<Pot>::duplicate(const Chem::Molecule& m) const
{
    bool res = false;
    for (std::size_t i = 0; i < conformers.size(); ++i) { // check geometry
        res = Numlib::kabsch_rmsd(conformers[i].xyz, m.get_xyz()) <= xtol;
        if (res) { // duplicate geometry; check energy
            double ediff = std::abs(conformers[i].energy - m.elec().energy());
            res = ediff <= etol;
        }
    }
    return res;
}

template <class Pot>
void Chem::Mcmm<Pot>::new_conformer()
{
    const unsigned ntrials = 20;
    // Generate a new random conformer by using the uniform usage scheme:
    for (unsigned i = 0; i < ntrials; ++i) {
        Numlib::Mat<double> xnew = mol.get_xyz();
        uniform_usage(xnew);
        mol.set_xyz(xnew);
        gen_rand_conformer(mol);
        if (accept_geom_dist(mol)) { // check geometry constraints
            break;
        }
    }

    // Perform geometry optimization:
    pot.run(mol);

    // Check acceptance:
    if (accept_energy(mol.elec().energy())) {
        if (!duplicate(mol)) { // store new conformer
            xcurr = mol.get_xyz();
            ecurr = mol.elec().energy();
            save_conformer(mol);
            naccept += 1;
        }
    }
}

template <class Pot>
void Chem::Mcmm<Pot>::update()
{
    kiter += 1;

    // Update global minimum if appropriate:
    if (!eglobal.empty()) {
        if (ecurr <= *std::min_element(eglobal.begin(), eglobal.end())) {
            eglobal.push_back(ecurr);
            xglobal = xcurr;
        }
    }
    else {
        eglobal.push_back(ecurr);
        xglobal = xcurr;
    }
    if (verbose) { // log state of MCMM solver
        double eglobal_min = *std::min_element(eglobal.begin(), eglobal.end());
        double ediff = std::abs(eglobal.end()[-1] - eglobal.end()[-2]);
        std::cout << "kiter = " << kiter << "; ecurr = " << ecurr
                  << "; eglobal = " << eglobal_min << "; conv = " << ediff
                  << std::endl;
    }
}

template <class Pot>
void Chem::Mcmm<Pot>::uniform_usage(Numlib::Mat<double>& xnew)
{
    xnew = xcurr;
    if (!conformers.empty()) {
        std::vector<Chem::Conformer> index;
        int istart = 0;
        int min_nstart = kiter;
        // Find least used:
        for (std::size_t i = 0; i < conformers.size(); ++i) {
            if (conformers[i].iter <= min_nstart) {
                index.push_back(conformers[i]);
            }
        }
        if (!conformers.empty()) { // find least used with lowest energy
            auto res = std::max_element(conformers.begin(), conformers.end());
            auto it = std::distance(conformers.begin(), res);
            double emin_ = conformers[it].energy;
            for (std::size_t i = 0; i < index.size(); ++i) {
                if (index[i].energy < emin_) {
                    emin_ = index[i].energy;
                    istart = i;
                }
            }
        }
        conformers[istart].iter += 1;
        xnew = conformers[istart].xyz;
    }
}

template <class Pot>
inline void Chem::Mcmm<Pot>::sort_conformers()
{
    std::sort(conformers.begin(), conformers.end());
    if (conformers.size() > nminima) {
        std::vector<Conformer> tmp;
        for (std::size_t i = 0; i < nminima; ++i) {
            tmp.push_back(Conformer(conformers[i]));
        }
        conformers.resize(nminima);
        conformers = tmp;
    }
}

template <class Pot>
std::vector<int> Chem::Mcmm<Pot>::select_rand_dihedral(const Chem::Molecule& m)
{
    auto connect = m.geom().get_connectivities();
    std::uniform_int_distribution<> rnd_uni_int(2, connect.size() - 1);
    int index = rnd_uni_int(mt);
    auto dihedral = connect[index];

    std::vector<int> res(0);
    for (std::size_t i = 2; i < connect.size(); ++i) {
        if (connect[i] == dihedral) {
            res.push_back(i);
        }
    }
    return res;
}

template class Chem::Mcmm<Chem::Gaussian>;
template class Chem::Mcmm<Chem::Mopac>;

