// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/molecule.h>
#include <numlib/constants.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_vibration")
{
    using namespace Chem;
    using namespace Numlib;
    using namespace Stdutils;

    SECTION("H2O")
    {
        std::ifstream from;
        fopen(from, "test_h2o.inp");

        Molecule mol(from);

        double zpe_ans = 0.024386;
        double zpe = mol.zero_point_energy() / Constants::au_to_icm;

        CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);

        Numlib::Vec<double> nu_ans = {2169.7566, 4141.6012, 4392.7809};
        Numlib::Vec<double> mu_ans = {1.0785, 1.0491, 1.0774};
        Numlib::Vec<double> k_ans = {2.9915, 10.6023, 12.2488};

        auto nu = mol.frequencies();
        auto mu = mol.vib_red_masses();
        auto k = mol.vib_force_constants();

        for (Index i = 0; i < nu.size(); ++i) {
            CHECK(std::abs(nu(i) - nu_ans(i)) < 1.0e-4);
        }
        for (Index i = 0; i < mu.size(); ++i) {
            CHECK(std::abs(mu(i) - mu_ans(i)) < 1.0e-4);
        }
        for (Index i = 0; i < k.size(); ++i) {
            CHECK(std::abs(k(i) - k_ans(i)) < 1.0e-4);
        }
    }

    SECTION("CO2")
    {
        std::ifstream from;
        fopen(from, "test_co2.inp");

        Molecule mol(from);

        double zpe_ans = 0.011626;
        double zpe = mol.zero_point_energy() / Constants::au_to_icm;
        CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);

        Numlib::Vec<double> nu_ans = {566.2442, 566.2442, 1435.2233, 2535.7089};
        Numlib::Vec<double> mu_ans = {12.8774, 12.8774, 15.9949, 12.8774};
        Numlib::Vec<double> k_ans = {2.4327, 2.4327, 19.4120, 48.7838};

        auto nu = mol.frequencies();
        auto mu = mol.vib_red_masses();
        auto k = mol.vib_force_constants();

        for (Index i = 0; i < nu.size(); ++i) {
            CHECK(std::abs(nu(i) - nu_ans(i)) < 1.0e-4);
        }
        for (Index i = 0; i < mu.size(); ++i) {
            CHECK(std::abs(mu(i) - mu_ans(i)) < 1.0e-4);
        }
        for (Index i = 0; i < k.size(); ++i) {
            CHECK(std::abs(k(i) - k_ans(i)) < 1.0e-4);
        }
    }

    SECTION("CH4OH")
    {
        std::ifstream from;
        fopen(from, "test_ch4oh.inp");

        Molecule mol(from);

        double zpe_ans = 0.051319;
        double zpe = mol.zero_point_energy() / Constants::au_to_icm;
        CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);

        Numlib::Vec<double> nu_ans = {
            -1752.6530, 62.301,    326.8547,  332.5309,  775.5512,
            909.5669,   1212.3977, 1275.3228, 1340.0812, 1451.6446,
            1476.6141,  3103.4053, 3244.1604, 3247.8661, 3767.9277};
        Numlib::Vec<double> mu_ans = {1.1201, 1.0445, 1.0905, 1.0684, 2.1719,
                                      1.5377, 1.1215, 1.0879, 1.1208, 1.0421,
                                      1.0216, 1.0242, 1.1068, 1.1074, 1.0669};
        Numlib::Vec<double> k_ans = {2.0272, 0.0024, 0.0686, 0.0696, 0.7697,
                                     0.7495, 0.9713, 1.0425, 1.1859, 1.2939,
                                     1.3124, 5.8116, 6.8635, 6.8826, 8.9245};

        auto nu = mol.frequencies();
        auto mu = mol.vib_red_masses();
        auto k = mol.vib_force_constants();

        for (Index i = 0; i < nu.size(); ++i) {
            CHECK(std::abs(nu(i) - nu_ans(i)) < 1.0e-4);
        }
        for (Index i = 0; i < mu.size(); ++i) {
            CHECK(std::abs(mu(i) - mu_ans(i)) < 1.0e-4);
        }
        for (Index i = 0; i < k.size(); ++i) {
            CHECK(std::abs(k(i) - k_ans(i)) < 1.0e-4);
        }
    }
}

