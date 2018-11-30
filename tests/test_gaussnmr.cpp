// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <fstream>
#include <string>
#include <cmath>

TEST_CASE("test_gaussnmr")
{
    std::ifstream from;
    Stdutils::fopen(from, "test_gaussnmr.inp");

    std::string nmr_method = "MP2 GIAO";
    double degen_tol = 0.05;

    std::vector<Chem::Gauss_NMR> nmr;
    std::vector<double> shield;

    Chem::Gauss_data gauss(from, Chem::out);
    gauss.get_nmr_data(nmr, nmr_method, degen_tol);

    for (auto& i : nmr) {
        shield.push_back(
            std::accumulate(i.shield.begin(), i.shield.end(), 0.0) /
            i.shield.size());
    }
    CHECK(shield.size() == 3);
    CHECK(std::abs(shield[0] - 31.6691) < 1.0e-12);
    CHECK(std::abs(shield[1] - 197.6891) < 1.0e-12);
    CHECK(std::abs(shield[2] - 362.2317) < 1.0e-12);
}

