// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <chem/collision.h>
#include <chem/thermodata.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program for computing the biased random walk model of R. G. Gilbert as
// published in the papers:
//
// "Collision energy exchange in highly vibrationally excited molecules:
// The biased random walk model". J. Chem. Phys., 1984, vol. 80, p. 5501.
//
// "Modeling collisional energy transfer in highly excited molecules"
// J. Chem. Phys., 1990, vol. 92, p. 1819.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Gilbert's biased random walk model");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>());
    // clang-format on

    auto result = options.parse(argc, argv);

    std::string input_file;

    if (result.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (result.count("file")) {
        input_file = result["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        std::ofstream to;

        std::string output_file;
        output_file = Stdutils::strip_suffix(input_file, ".inp");
        output_file = output_file + ".out";

        Stdutils::fopen(from, input_file);
        Stdutils::fopen(to, output_file.c_str());

        Chem::Thermodata tpdata(from);
        auto temp = tpdata.get_temperature();

        Chem::Collision model(from);
        for (auto& t : temp) {
            model.biased_random_walk(t, to);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

