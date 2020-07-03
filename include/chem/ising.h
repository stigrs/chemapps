// Copyright (c) 2020 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_ISING_H
#define CHEM_ISING_H

#include <numlib/matrix.h>
#include <numlib/math.h>
#include <array>
#include <random>

namespace Chem {

// Class providing the two-dimensional Ising model.
//
class Ising2D {
public:
    Ising2D(int sz, double j, double b, int seed = 0);

    // Perform Metropolis algorithm.
    std::array<double, 4> metropolis(double temp, int mc_trials = 1000);

private:
    // Initialise spins for the ground state.
    void init_spins();

    // Perform Monte Carlo sampling of spin flips.
    void mc_spin_flip(double beta);

    // Compute lattice energy.
    void compute_energy_magn();

    // Compute magnetisation.
    double magnetisation() const;

    // Periodic boundary conditions.
    int pbc(int i) const;

    int size;      // lattice size
    double jint;   // interaction (ferromagnetic if positive)
    double bfield; // external magnetic field

    double energy; // lattice energy
    double magn;   // net magnetisation

    Numlib::Mat<int> spins; // spin matrix

    std::mt19937_64 mt; // random number engine
};

inline Ising2D::Ising2D(int sz, double j, double b, int seed)
    : size(sz), jint(j), bfield(b)
{
    if (seed == 0) {
        std::random_device rd;
        std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(ss);
    }
    else {
        mt.seed(seed); // should only be used for testing
    }
    init_spins();
}

inline void Ising2D::init_spins()
{
    spins = Numlib::ones<Numlib::Mat<int>>(size, size);
}

inline double Ising2D::magnetisation() const
{
    return static_cast<double>(Numlib::sum(spins.ravel()));
}

inline int Ising2D::pbc(int i) const { return (i + size) % size; }
} // namespace Chem

#endif // CHEM_ISING_H
