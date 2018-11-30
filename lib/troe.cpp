// Copyright (c) 2013-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/troe.h>
#include <chem/whirab.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <stdexcept>
#include <cmath>

Chem::Troe::Troe(std::istream& from, Molecule& mol_) : mol(mol_)
{
    using namespace Stdutils;

    int pot_type_tmp = 1;
    auto pos = find_token(from, "Troe");
    if (pos != -1) {
        get_token_value(from, pos, "pot_type", pot_type_tmp, pot_type_tmp);
        get_token_value(from, pos, "en_barrier", en_barrier, 0.0);
        get_token_value(from, pos, "imom_ratio", imom_ratio, 1.0);
        get_token_value(from, pos, "n_free_rot", n_free_rot, 0);
        get_token_value(from, pos, "n_morse_osc", n_morse_osc, 0);
    }
    // Check if data are sensible:

    if (pot_type_tmp == 1) {
        pot_type = type1;
    }
    else if (pot_type_tmp == 2) {
        pot_type = type2;
    }
    else {
        throw std::runtime_error("bad potential type");
    }
    Assert::dynamic(en_barrier > 0.0, "bad energy barrier");
    Assert::dynamic(imom_ratio > 0.0, "bad moment of inertia ratio");
    Assert::dynamic(n_free_rot >= 0, "bad number of free rotors");
    Assert::dynamic(n_morse_osc >= 0, "bad number of Morse oscillators");

    zpe = mol.vib().zero_point_energy();
}

double Chem::Troe::f_energy(const double temp) const
{
    using namespace Numlib::Constants;

    double en = en_barrier + Whirab::a_corr(mol, en_barrier) * zpe;
    en = R * 1.0e-3 * temp / (en * icm_to_kJ);

    auto s = mol.vib().frequencies().size();
    double f_e = 0.0;

    for (Index i = 0; i < s; ++i) {
        f_e += std::tgamma(s) * std::pow(en, i) / std::tgamma(s - i);
    }
    return f_e;
}

double Chem::Troe::f_rotation(const double temp) const
{
    using namespace Numlib::Constants;
    const double kT = R * 1.0e-3 * temp / icm_to_kJ;

    double e0_azpe = en_barrier + Whirab::a_corr(mol, en_barrier) * zpe;
    double e0_kT = en_barrier / kT;
    double f_rot = 1.0;
    auto s = mol.vib().frequencies().size();

    switch (pot_type) {
    case type1: // no barrier
        if (mol.structure() == linear) {
            e0_azpe /= s * kT;
            e0_kT = 2.15 * std::pow(e0_kT, 1. / 3.);
            f_rot = e0_azpe * (e0_kT / (e0_kT - 1.0 + e0_azpe));
        }
        else { // nonlinear
            e0_kT = std::pow(e0_kT, 1. / 3.);
            f_rot = (std::tgamma(s) / std::tgamma(s + 0.5 + 1.0)) *
                    std::pow(e0_azpe / kT, 1.5) *
                    (2.15 * e0_kT /
                     (2.15 * e0_kT - 1.0 + e0_azpe / ((s + 0.5) * kT)));
        }
        break;
    case type2: // barrier
        if (mol.structure() == linear) {
            e0_azpe /= s * kT;
            f_rot = e0_azpe * imom_ratio / (imom_ratio - 1.0 - e0_azpe);
        }
        else { // nonlinear
            f_rot =
                (std::tgamma(s) / std::tgamma(s + 0.5 + 1.0)) *
                std::pow(e0_azpe / kT, 1.5) *
                (imom_ratio / (imom_ratio - 1.0 + e0_azpe / ((s + 0.5) * kT)));
        }
        break;
    }
    return f_rot;
}

double Chem::Troe::f_free_rotor(const double temp) const
{
    using namespace Numlib::Constants;
    double f_free_rot = 1.0;

    if (n_free_rot > 0) {
        double kT = R * 1.0e-3 * temp / icm_to_kJ;
        double e0_azpe =
            (en_barrier + Whirab::a_corr(mol, en_barrier) * zpe) / kT;

        auto s = mol.vib().frequencies().size();
        auto r = n_free_rot;

        f_free_rot = (std::tgamma(s) / std::tgamma(s + 0.5 * r)) *
                     std::pow(e0_azpe, 0.5 * r);
    }
    return f_free_rot;
}

double Chem::Troe::f_hind_rotor(const double temp) const
{
    using namespace Numlib::Constants;
    double f_hind_rot = 1.0;

    if (mol.tor().pot_coeff().size() > 0) {
        const double kT = k * temp;
        const double f = icm_to_kJ * 1.0e+3 / N_A;

        double v0 = Numlib::max(mol.tor().pot_coeff());
        if ((en_barrier / v0) <= 3.0) {
            throw std::runtime_error(
                "Troe::f_hind_rotor(): E0/V0 <= 3, not implemented yet");
        }
        double a = Whirab::a_corr(mol, en_barrier);
        double s = mol.vib().frequencies().size();
        double n = mol.tor().symmetry_number();
        double b = Numlib::max(mol.tor().constant()) * 100.0;
        double v0f = v0 * f;
        double e0 = en_barrier * f;
        double zpef = zpe * f;

        double denom =
            std::pow((1.0 - exp(-kT / v0f)), 1.2) +
            exp(-1.2 * kT / v0f) /
                (std::sqrt(kT / (2.0 * b * c_0 * h_bar)) *
                 (1.0 -
                  std::exp(-std::sqrt(n * n * h * c_0 * b * v0f / (kT * kT)))));

        f_hind_rot = (std::tgamma(s) / std::tgamma(s - 0.5 + 1.0));
        f_hind_rot *= std::sqrt((e0 + a * zpef) / kT);
        f_hind_rot *= 1.0 - std::exp(-e0 / (s * v0f));
        f_hind_rot /= denom;
    }
    return f_hind_rot;
}

