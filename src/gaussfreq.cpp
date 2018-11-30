// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

//  Extracts vibrational frequencies from Gaussian output file.
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

    std::vector<double> freqs;

    Chem::Gauss_data gauss(from, Chem::out);
    gauss.get_freqs(freqs);

    for (auto vi : freqs) {
        std::cout << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4) << vi << '\t';
    }
    std::cout << '\n';
}
