// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>

//  Extracts optimized energies and geometry from Gaussian output file.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.log\n";
        return 1;
    }

    std::ifstream from(args[1]);
    if (!from) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }

    Chem::Gauss_data gauss(from, Chem::out);
    auto en = gauss.get_scf_zpe_energy();
    double en_sum = std::accumulate(en.begin(), en.end(), 0.0);

    Chem::Gauss_coord coord;
    gauss.get_opt_cart_coord(coord);

    Stdutils::Format<double> fix;
    fix.fixed().width(15).precision(8);
    std::cout << "SCF: " << fix(en[0]) << " Hartree\n"
              << "ZPE: " << fix(en[1]) << " Hartree\n"
              << "Tot: " << fix(en_sum) << " Hartree\n\n";

    gauss.print_opt_geom();
}
