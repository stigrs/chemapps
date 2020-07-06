// Copyright (c) 2020 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <chem/ising.h>
#include <cxxopts.hpp>
#include <exception>
#include <iostream>
#include <cmath>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program providing two-dimensional Ising solver.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Two-dimensional Ising solver");
    options.add_options()
        ("h,help", "display help message")
        ("size", "lattice size", cxxopts::value<int>()) 
        ("jint", "interaction", cxxopts::value<double>()->default_value("1.0"))
        ("bfield", "external magnetic field", cxxopts::value<double>()->default_value("0.0")) 
        ("t0", "intial temperature", cxxopts::value<double>()) 
        ("t1", "final temperature", cxxopts::value<double>()) 
        ("ntemp", "number of temperature steps", cxxopts::value<int>()) 
        ("trials", "number of Monte Carlo trials", cxxopts::value<int>()->default_value("1000")) 
        ("viz", "vizualisation steps", cxxopts::value<std::vector<int>>());
    // clang-format on

    auto args = options.parse(argc, argv);

    int size = 0;
    int ntemp = 0;
    int trials = 0;

    double t0 = 0.0;
    double t1 = 0.0;
    double jint = 0.0;
    double bfield = 0.0;

    std::vector<int> viz;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("size")) {
        size = args["size"].as<int>();
    }
    if (args.count("ntemp")) {
        ntemp = args["ntemp"].as<int>();
    }
    if (args.count("t0")) {
        t0 = args["t0"].as<double>();
    }
    if (args.count("t1")) {
        t1 = args["t1"].as<double>();
    }
    if (args.count("viz")) {
        viz = args["viz"].as<std::vector<int>>();
    }
    trials = args["trials"].as<int>();
    jint = args["jint"].as<double>();
    bfield = args["bfield"].as<double>();

    try {
        Chem::Ising2D mod(size, jint, bfield, viz);
        std::cout << "T,E,<M>,Cv,X\n";

        auto temp = Numlib::linspace(t0, t1, ntemp);
        for (auto& ti : temp) {
            auto res = mod.metropolis(ti, trials);
            std::cout << ti << "," << res[0] << "," << std::abs(res[1]) << ","
                      << res[2] << "," << res[3] << std::endl;
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
