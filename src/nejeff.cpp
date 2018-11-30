// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <fstream>
#include <iostream>

// Program for calculation of effective N(E,J) = N1 * N2 / (N1 + N2).
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 3) {
        std::cerr << "usage: " << args[0] << " flux_file_1 flux_file_2\n";
        return 1;
    }

    std::ifstream from1(args[1]);
    if (!from1) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }
    std::ifstream from2(args[2]);
    if (!from2) {
        std::cerr << "cannot open " << args[2] << '\n';
        return 1;
    }

    double e1;
    double e2;
    double j1;
    double j2;
    double n1;
    double n2;
    double neff;

    while (from1 >> e1 >> j1 >> n1) {
        from2 >> e2 >> j2 >> n2;
        if (!from2) {
            std::cerr << "input flux file " << args[2] << " too short?\n";
            return 1;
        }
        if (e1 != e2) {
            std::cerr << "bad E grid: " << e1 << ", " << e2 << '\n';
            return 1;
        }
        if (j1 != j2) {
            std::cerr << "bad J grid: " << j1 << ", " << j2 << '\n';
            return 1;
        }
        if ((n1 <= 0.0) || (n2 <= 0.0)) {
            neff = 0.0;
        }
        else {
            neff = n1 * n2 / (n1 + n2);
        }
        std::cout << e1 << " " << j1 << " " << neff << '\n';
    }
    if (from2 >> n2) {
        std::cerr << "input flux file " << args[1] << " too short?\n";
        return 1;
    }
}

