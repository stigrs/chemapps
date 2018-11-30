// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/energy_levels.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <cmath>

std::vector<double> Chem::Energy_levels::harmonic_oscillator(double freq,
                                                             double emax)
{
    Assert::dynamic<Assert::level(2)>(freq > 0.0, "bad freq");
    Assert::dynamic<Assert::level(2)>(emax > 0.0, "bad emax");

    int kmax = 1 + Numlib::round<int>(emax / freq);
    std::vector<double> result(kmax);
    for (int i = 0; i < kmax; ++i) {
        result[i] = freq * (i + 1);
    }
    return result;
}

std::vector<double> Chem::Energy_levels::free_rotor(double rotc, double emax)
{
    Assert::dynamic<Assert::level(2)>(rotc > 0.0, "bad rotc");
    Assert::dynamic<Assert::level(2)>(emax > 0.0, "bad emax");

    std::vector<double> result;
    double ej = 0.0;
    int j = 1;
    while (ej < emax) {
        ej = rotc * j * j;
        result.push_back(ej);
        ++j;
    }
    return result;
}

std::vector<double> Chem::Energy_levels::hindered_rotor(double sigma,
                                                        double rotc,
                                                        double barrier,
                                                        double emax)
{
    Assert::dynamic<Assert::level(2)>(sigma >= 1.0, "bad sigma");
    Assert::dynamic<Assert::level(2)>(rotc > 0.0, "bad rotc");
    Assert::dynamic<Assert::level(2)>(barrier >= 0.0, "bad barrier");
    Assert::dynamic<Assert::level(2)>(emax > 0.0, "bad emax");

    std::vector<double> result;

    if (barrier > 1.0) { // hindered rotor
        double sig = sigma;
        double b = rotc;
        double v0 = barrier;
        double frq = sig * std::sqrt(b * v0);
        double r = v0 / frq;
        double es = 0.0;
        double zpe = 0.0;
        int ns = 0;
        while (es < emax) {
            // Calculate harmonic oscillator energy level:
            double dnv = ns / sig;
            int nv = Numlib::round<int>(dnv);
            double tv = -frq * (1.0 + 2.0 * nv + 2.0 * nv * nv) / (16.0 * r);
            double ev = frq * (nv + 0.5) + tv;

            // Calculate free rotor energy level:
            int j = (ns + 1) / 2;
            double tr = 0.0;
            if ((j > (r * sig / 2.0)) && (r > 0.0)) {
                tr = std::pow(r, 4.0) * sig * sig * b /
                     (8.0 * (std::pow(2.0 * j / sig, 2.0) - 1.0));
            }
            double ej = b * j * j + 0.5 * v0 + tr;

            // Calculate hindered rotor energy level:
            double s = 0.5 * (1.0 + std::tanh(5.0 * (ev - v0) / v0));
            if (ej > 1.5 * v0) {
                s = 1.0;
            }
            es = ev * (1.0 - s) + ej * s;

            if (ns == 0) {
                zpe = es;
            }
            if (ns > 0) {
                result.push_back(es - zpe);
            }
            ++ns;
        }
    }
    else { // free rotor
        result = free_rotor(rotc, emax);
    }
    return result;
}

