// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/thermodata.h>
#include <stdutils/stdutils.h>

Chem::Thermodata::Thermodata(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    // Read input data:

    Numlib::Vec<double> p_def = {Numlib::Constants::std_atm};
    Numlib::Vec<double> t_def = {298.15};

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "pressure", pressure, p_def);
        get_token_value(from, pos, "temperature", temperature, t_def);
        get_token_value(from, pos, "incl_sigma", incl_sigma, 1);
        get_token_value(from, pos, "zeroref", zeroref, std::string("BOT"));
    }

    // Validate input:

    for (const auto& p : pressure) {
        Assert::dynamic(p > 0.0, "bad pressure");
    }
    for (const auto& t : temperature) {
        Assert::dynamic(t > 0.0, "bad temperature");
    }
    Assert::dynamic(incl_sigma == 0 || incl_sigma == 1, "bad incl_sigma");
    Assert::dynamic(zeroref == "BOT" || zeroref == "V=0", "bad zeroref");
}

