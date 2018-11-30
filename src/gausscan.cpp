// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Program for extracting data from a Gaussian relaxed PES scan.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.log\n";
        return 1;
    }

    try {
        std::ifstream from;
        Stdutils::fopen(from, args[1]);

        std::string scan_coord;
        std::vector<double> coord;
        std::vector<double> energy;

        Chem::Gauss_data gauss(from, Chem::out);
        gauss.get_pes_scan_data(scan_coord, coord, energy);

        Stdutils::Format<double> fix6;
        fix6.fixed();

        std::cout << "#\t" << scan_coord << "\t\t"
                  << "Eigenvalues\n";
        for (std::size_t i = 0; i < energy.size(); ++i) {
            std::cout << i + 1 << '\t' << fix6(coord[i]) << '\t'
                      << fix6(energy[i]) << '\n';
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

