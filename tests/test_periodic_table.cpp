// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/periodic_table.h>
#include <catch2/catch.hpp>

TEST_CASE("test_periodic_table")
{
    using namespace Chem::Periodic_table;

    SECTION("Carbon")
    {
        double atomic_weight = get_element("C").atomic_weight;
        CHECK(atomic_weight == (12.0096 + 12.0116) * 0.5);
    }

    SECTION("109Ag")
    {
        double atomic_mass = get_element("109Ag").atomic_mass;
        CHECK(atomic_mass == 108.904755);
    }

    SECTION("atomic_number") { CHECK(get_atomic_symbol(15) == "P"); }
}
