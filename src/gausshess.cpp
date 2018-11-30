// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

//  Extracts Cartesian force constants from Gaussian fchk file.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.fchk\n";
        return 1;
    }
    std::ifstream from(args[1]);
    if (!from) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }

    const std::string pattern = "Cartesian Force Constants";

    bool found = false;
    double fc;
    std::string line;
    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            found = true;
            int count = 0;
            std::cout << " HESSIAN\n";
            while (from >> fc) {
                std::cout << std::setw(16) << std::setprecision(8)
                          << std::setiosflags(std::ios_base::scientific |
                                              std::ios_base::uppercase)
                          << fc;
                count++;
                if (count == 5) {
                    std::cout << '\n';
                    count = 0;
                }
            }
            if (count != 0) {
                std::cout << '\n';
            }
            std::cout << " END\n\n";
        }
    }
    if (!found) {
        std::cerr << "could not find force constants\n";
        return 1;
    }
}

