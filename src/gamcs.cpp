// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#pragma warning(disable : 4996)      // caused by ctime
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
#include <chrono>

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

        Stdutils::fopen(from, input_file);
        Stdutils::fopen(to, output_file);

        auto tstart = std::chrono::system_clock::now();

        if (pot == "Gaussian" || pot == "gaussian") {
            Chem::Gamcs<Chem::Gaussian> ga(from, to);
            ga.solve(to);
        }
        else {
            Chem::Gamcs<Chem::Mopac> ga(from, to);
            ga.solve(to);
        }

        auto tend = std::chrono::system_clock::now();
        std::chrono::duration<double> telapsed = tend - tstart;
        std::time_t end_time = std::chrono::system_clock::to_time_t(tend);

        to << "Computation finished at " << std::ctime(&end_time) << '\n'
           << "Elapsed time: " << telapsed.count() << " s\n";
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif
