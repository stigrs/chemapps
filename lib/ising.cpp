// Copyright (c) 2020 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/ising.h>
#include <stdexcept>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

std::array<double, 4> Chem::Ising2D::metropolis(double temp, int mc_trials)
{
    // Initialise spin matrix and variables:
    auto spins = init_spins();

    double beta = 1.0 / temp;
    double e1 = 0.0;
    double e2 = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;

    // Perform equilibration:
    for (int it = 0; it < mc_trials; ++it) {
        mcmove(spins, beta);
    }

    // Perform Monte Carlo sampling:
    for (int it = 0; it < mc_trials; ++it) {
        mcmove(spins, beta);
        compute_energy_magn(spins);
        e1 += energy;
        m1 += magn;
        e2 += energy * energy;
        m2 += magn * magn;
    }
    // Normalise average values:
    e1 /= static_cast<double>(mc_trials);
    e2 /= static_cast<double>(mc_trials);
    m1 /= static_cast<double>(mc_trials);
    m2 /= static_cast<double>(mc_trials);
    double n2 = static_cast<double>(size * size);

    double e_avg = e1 / n2;
    double m_avg = m1 / n2;
    double c_avg = beta * beta * (e2 - e1 * e1) / n2;
    double x_avg = beta * (m2 - m1 * m1) / n2;

    return {e_avg, m_avg, c_avg, x_avg};
}

void Chem::Ising2D::mcmove(Numlib::Mat<int>& spins, double beta)
{
    std::uniform_int_distribution<> rnd_int_uni(0, size - 1);
    std::uniform_real_distribution<> rnd_real_uni(0.0, 1.0);

    for (int it = 0; it < size * size; ++it) {
        int i = rnd_int_uni(mt);
        int j = rnd_int_uni(mt);

        auto st = spins(i, j);
        auto s_nb = spins(i, pbc(j - 1)) + spins(i, pbc(j + 1)) +
                    spins(pbc(i - 1), j) + spins(pbc(i + 1), j);

        // Compute energy change:
        double ediff = 2.0 * bfield * st + 2.0 * jint * st * s_nb;

        // Check for acceptance; flip spin if accepted:
        if (ediff < 0.0) {
            st *= -1;
        }
        else if (rnd_real_uni(mt) <= std::exp(-ediff * beta)) {
            st *= -1;
        }
        spins(i, j) = st;
    }
}

void Chem::Ising2D::compute_energy_magn(const Numlib::Mat<int>& spins)
{
    magn = magnetisation(spins);
    energy = 0.0;
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            auto sij = spins(i, j);
            auto s_nb = spins(i, pbc(j + 1)) + spins(pbc(i + 1), j);
            energy -= sij * s_nb;
        }
    }
    energy *= jint;
    energy -= bfield * magn;
}
