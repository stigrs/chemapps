// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <chem/gaussian.h>
#include <chem/mopac.h>
#include <chem/gamcs.h>
#include <chem/molecule.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program providing Genetic Algorithm Molecular Structure Search.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Genetic Algorithm Molecular Conformer Search");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>()) 
		("p,pot", "potential (Gaussian or Mopac)", cxxopts::value<std::string>());
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;
    std::string pot = "Mopac";

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("pot")) {
        pot = args["pot"].as<std::string>();
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        std::ofstream to;

        std::string output_file =
            Stdutils::strip_suffix(input_file, "inp") + "out";
        std::cout << output_file << std::endl;

        Stdutils::fopen(from, input_file);
        Stdutils::fopen(to, output_file);

        if (pot == "Gaussian" || pot == "gaussian") {
            Chem::Gamcs<Chem::Gaussian> ga(from, to);
            ga.solve(to);
        }
        else {
            Chem::Gamcs<Chem::Mopac> ga(from, to);
            ga.solve(to);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

