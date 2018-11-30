// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/thermodata.h>
#include <chem/tunnel.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <fstream>

TEST_CASE("test_tunnel")
{
    SECTION("eckart")
    {
        std::ifstream from;
        Stdutils::fopen(from, "test_tunnel.inp");

        Chem::Tunnel tunnel(from);
        Chem::Thermodata td(from);

        // Results from Brown (1981):
        Numlib::Vec<double> ans = {6.41025, 2.42735, 1.29311, 1.17025, 1.11455};
        auto temp = td.get_temperature();

        for (int i = 0; i < temp.size(); ++i) {
            double res = tunnel.factor(temp(i));
            CHECK(std::abs(res - ans(i)) < 1.0e-3);
        }
    }

    SECTION("none")
    {
        Chem::Tunnel tunnel_none;
        CHECK(tunnel_none.factor() == 1.0);
    }
}

