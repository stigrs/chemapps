// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/collision.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <fstream>
#include <cmath>

TEST_CASE("test_collision")
{
    SECTION("C3H7Br-Ne")
    {
        // Data taken from UNIMOL program of Gilbert.

        std::ifstream from;
        Stdutils::fopen(from, "test_c3h7br_ne.inp");

        Chem::Collision coll(from);
        CHECK(std::abs(coll.sigma_complex() - 4.30) < 1.0e-12);
        CHECK(std::abs(coll.epsilon_complex() - 113.999) < 1.0e-3);
        CHECK(std::abs(coll.sigma_local() - 2.96) < 1.0e-2);
        CHECK(std::abs(coll.epsilon_local() - 29.36) < 1.0e-2);
        CHECK(std::abs(coll.reduced_mass() - 17.183) < 1.0e-3);
        CHECK(std::abs(coll.average_mass() - 10.080) < 1.0e-3);
        CHECK(std::abs(coll.impact_parameter(870.0) - 2.867) < 1.0e-3);
        CHECK(std::abs(coll.collision_time(870.0) * 1.0e+13 - 1.45) < 1.0e-3);
        CHECK(std::abs(coll.s_parameter(870.0) - 542.986) < 1.0e-1);
    }
}
