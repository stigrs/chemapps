// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/whirab.h>
#include <chem/traits.h>
#include <numlib/constants.h>
#include <numlib/traits.h>
#include <cmath>

double Chem::Whirab::a_corr(const Chem::Molecule& mol, double e_barrier)
{
    double en = e_barrier / mol.vib().zero_point_energy();
    double w;
    if (en >= 1.0) {
        w = -1.0506 * std::pow(en, 0.25);
        w = std::pow(10.0, w);
    }
    else {
        w = 1.0 / (5.0 * en + 2.73 * std::sqrt(en) + 3.51);
    }
    double sum2_v = 0.0;
    double sum_v2 = 0.0;
    for (auto vi : mol.vib().frequencies()) {
        if (vi < 0.0) { // ignore imaginary frequencies
            continue;
        }
        else {
            sum2_v += vi;
            sum_v2 += vi * vi;
        }
    }
    sum2_v *= sum2_v;

    double factor = 3.0; // structure factor for nonlinear
    if (mol.structure() == linear) {
        factor = 2.0;
    }
    auto s = narrow_cast<double>(mol.vib().frequencies().size());
    auto r = narrow_cast<double>(mol.tor().pot_coeff().size());

    double beta =
        (s - 1.0) * ((s + 0.5 * r + 0.5 * factor) / s) * sum_v2 / sum2_v;

    return 1.0 - beta * w;
}

double Chem::Whirab::vibr_density_states(const Chem::Molecule& mol,
                                         double e_barrier)
{
    double hv = 1.0;
    for (auto wi : mol.vib().frequencies()) {
        if (wi < 0.0) { // ignore imaginary frequencies
            continue;
        }
        else {
            hv *= wi;
        }
    }
    double s = narrow_cast<double>(mol.vib().frequencies().size());
    double sm1 = s - 1.0;
    double rho =
        e_barrier + a_corr(mol, e_barrier) * mol.vib().zero_point_energy();

    rho = std::pow(rho, sm1) / (std::tgamma(s) * hv);
    return rho / Numlib::Constants::icm_to_kJ;
}

