// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program for generating XYZ files from MCMM solver output.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Generate XYZ files from MCMM solver output");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>()) 
		("N,atoms", "number of atoms", cxxopts::value<int>()) 
        ("t,title", "title line", cxxopts::value<std::string>());
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;
    std::string title = "Title";
    int natoms = 0;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("title")) {
        title = args["title"].as<std::string>();
    }
    if (args.count("atoms")) {
        natoms = args["atoms"].as<int>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
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
        std::string line;
        while (std::getline(from, line)) {
            if (line.find("Conformer:") != std::string::npos) {
                std::string token;
                int nc;
                std::istringstream iss(line);
                iss >> token >> nc;

                std::ofstream to;
                std::string base = Stdutils::strip_suffix(input_file, ".out");
                base += "_c" + std::to_string(nc);
                std::string com = base + ".xyz";
                Stdutils::fopen(to, com.c_str());

                to << natoms << '\n' << title << '\n';

                for (int it = 0; it < 5; ++it) {
                    std::getline(from, line); // ignore
                }
                int center;
                std::string atom;
                double x;
                double y;
                double z;
                while (std::getline(from, line)) {
                    iss = std::istringstream(line);
                    if (line.find("--") != std::string::npos) {
                        break;
                    }
                    if (iss >> center >> atom >> x >> y >> z) {
                        Stdutils::Format<double> fix;
                        fix.fixed().width(10).precision(6);
                        to << atom << '\t' << fix(x) << "  " << fix(y) << "  "
                           << fix(z) << '\n';
                    }
                }
                to << '\n';
            }
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

