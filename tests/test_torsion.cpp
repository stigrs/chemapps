// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/molecule.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_torsion")
{
    using namespace Chem;
    using namespace Stdutils;

    SECTION("CH2ClCH2Cl")
    {
        std::ifstream from;
        fopen(from, "test_ch2clch2cl.inp");

        Molecule mol(from);
        double rmi = mol.tor().red_moment();

        const double rmi_ans = 58.76991427; // Chuang and Truhlar (2000)
        CHECK(std::abs(rmi - rmi_ans) < 9.0e-2);
    }

    SECTION("CH3OH")
    {
        std::ifstream from;
        fopen(from, "test_ch3oh.inp");

        Molecule mol(from);
        double rmi = mol.tor().red_moment();

        const double rmi_ans = 2.19; // Chuang and Truhlar (2000)
        CHECK(std::abs(rmi - rmi_ans) < 2.0e-2);
    }
}
