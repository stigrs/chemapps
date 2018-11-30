// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/electronic.h>
#include <stdutils/stdutils.h>

Chem::Electronic::Electronic(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    charge_ = 0;
    spin = 1;
    energy_ = 0.0;
    so_degen = {1};
    so_energy = {0.0};

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "charge", charge_, 0);
        get_token_value(from, pos, "spin_mult", spin, 1);
        get_token_value(from, pos, "elec_energy", energy_, 0.0);
        get_token_value(from, pos, "so_degen", so_degen, so_degen);
        get_token_value(from, pos, "so_energy", so_energy, so_energy);
    }
    Assert::dynamic(same_extents(so_degen, so_energy), "bad spin-orbit input");
}

