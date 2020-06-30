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
    Numlib::Mat<int> init_spins() const;

    // Perform Monte Carlo sampling.
    void mcmove(Numlib::Mat<int>& spins, double beta);

    // Compute lattice energy.
    double energy(const Numlib::Mat<int>& spins) const;

    // Compute magnetisation.
    double magnetisation(const Numlib::Mat<int>& spins) const;

	// Periodic boundary conditions.
    int pbc(int i) const;

    int size;      // lattice size
    double jint;   // interaction (ferromagnetic if positive)
    double bfield; // external magnetic field

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
}

inline Numlib::Mat<int> Ising2D::init_spins() const
{
    return Numlib::ones<Numlib::Mat<int>>(size, size);
}

inline double Ising2D::magnetisation(const Numlib::Mat<int>& spins) const
{
    return static_cast<double>(Numlib::sum(spins.ravel()));
}

inline int Ising2D::pbc(int i) const
{
    return (i + size) % size;
}
} // namespace Chem

#endif // CHEM_ISING_H
