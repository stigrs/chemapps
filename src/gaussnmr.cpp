// Copyright (c) 2013-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

//  Summarize NMR shieldings from a Gaussian calculation.
//  The program produces the same output as GaussView.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Summarize Gaussian NMR calculation");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>())
        ("m,method", "NMR method (default is GIAO)", cxxopts::value<std::string>())
        ("t,tol", "degeneracy tolerance (default is 0.05)", cxxopts::value<double>());
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;
    std::string nmr_method = "SCF GIAO"; // default NMR method
    double degen_tol = 0.05;             // default degeneracy tolerance

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }
    if (args.count("method")) {
        nmr_method = args["method"].as<std::string>();
    }
    if (args.count("tol")) {
        degen_tol = args["tol"].as<double>();
    }

    std::ifstream from;
    Stdutils::fopen(from, input_file);

    // Get NMR data:

    std::vector<Chem::Gauss_NMR> nmr;

    Chem::Gauss_data gauss(from, Chem::out);
    gauss.get_nmr_data(nmr, nmr_method, degen_tol);

    // Output results:

    Stdutils::Format<char> hline;
    hline.width(37).fill('-');

    std::cout << "\nSummary of NMR spectrum (" << nmr_method
              << " magnetic shieldings)\n"
              << "Degenerate peaks are condensed together "
              << "(degeneracy tolerance " << degen_tol << ")\n\n"
              << "Shielding/ppm\t"
              << "Degen.\t"
              << "Elem.\t"
              << "Atoms\n"
              << hline('-') << '\n';

    Stdutils::Format<double> fix8;
    fix8.fixed().width(9).precision(4);

    for (auto& ni : nmr) {
        double shield =
            std::accumulate(ni.shield.begin(), ni.shield.end(), 0.0) /
            static_cast<double>(ni.shield.size());

        std::cout << fix8(shield) << '\t' << ni.shield.size() << '\t' << ni.atom
                  << '\t';

        std::sort(ni.number.begin(), ni.number.end());

        for (std::size_t j = 0; j < ni.number.size(); ++j) {
            std::cout << ni.number[j];
            if ((ni.number.size() > 1) && (j < ni.number.size() - 1)) {
                std::cout << ',';
            }
        }
        std::cout << '\n';
    }
}

