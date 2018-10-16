// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/impl/elec_state.h>
#include <stdutils/stdutils.h>
#include <cassert>

Chem::Impl::Elec_state::Elec_state(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    so_degen = {1};
    so_energy = {0.0};

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "charge", charge, 0);
        get_token_value(from, pos, "spin_mult", spin, 1);
        get_token_value(from, pos, "elec_energy", energy, 0.0);
        get_token_value(from, pos, "so_degen", so_degen, so_degen);
        get_token_value(from, pos, "so_energy", so_energy, so_energy);
    }
    assert(same_extents(so_degen, so_energy));
}

