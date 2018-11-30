// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/units.h>
#include <catch2/catch.hpp>

TEST_CASE("units")
{
    Chem::Units::Type ans = Chem::Units::Type::kJ_mol;
    CHECK(Chem::Units::lexer("kJ/mol") == ans);
}

