// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/molecule.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_rotation")
{
    using namespace Chem;
    using namespace Numlib;
    using namespace Stdutils;

    std::ifstream from;
    fopen(from, "test_rotation.inp");
    Molecule mol(from);

    Numlib::Vec<double> ans = {127.63201, 24.89071, 24.02767};
    Numlib::Vec<double> res = mol.rot().constants();

    for (Index i = 0; i < ans.size(); ++i) {
        CHECK(std::abs(res(i) - ans(i)) < 1.0e-4);
    }
}
