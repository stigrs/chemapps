// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/thermodata.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <string>
#include <cmath>

TEST_CASE("test_thermodata")
{
    Numlib::Vec<double> t_ans = {100.0, 200.0, 300.0, 800.0, 1000.0};
    std::string zeroref_ans = "V=0";

    std::ifstream from;
    Stdutils::fopen(from, "test_thermodata.inp");

    Chem::Thermodata td(from);

    CHECK(td.get_vibr_zeroref() == zeroref_ans);
    for (Index i = 0; i < t_ans.size(); ++i) {
        CHECK(std::abs(td.get_temperature()(i) - t_ans(i)) < 1.0e-12);
    }
}

