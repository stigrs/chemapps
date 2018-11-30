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
#include <chem/mcmm.h>
#include <chem/molecule.h>
#include <chem/mopac.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program providing Monte Carlo Multiple Minima (MCMM) solver.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Monte Carlo Multiple Minima (MCMM) solver");
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
        Stdutils::fopen(from, input_file);

        Chem::Molecule mol(from);
        if (pot == "Gaussian" || pot == "gaussian") {
            Chem::Mcmm<Chem::Gaussian> mc(from, mol, "Mcmm", true);
            mc.solve();
        }
        else {
            Chem::Mcmm<Chem::Mopac> mc(from, mol, "Mcmm", true);
            mc.solve();
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

