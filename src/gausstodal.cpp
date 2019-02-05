// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gauss_data.h>
#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>
#include <iostream>
#include <fstream>
#include <vector>

// Get Cartesian coordinates from Gaussian output file and convert to
// Dalton mol file format.
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
    // Get optimized Cartesian coordinates:
    Chem::Gauss_data gauss(from, Chem::out);
    Chem::Gauss_coord data;
    gauss.get_opt_cart_coord(data);

    // Count number of different atom types:
    std::vector<int> atoms;
    atoms.push_back(data.atnum[0]);

    for (std::size_t i = 1; i < data.atnum.size(); ++i) {
        if (std::find(atoms.begin(), atoms.end(), data.atnum[i]) ==
            atoms.end()) {
            atoms.push_back(data.atnum[i]);
        }
    }
    // Count number of atoms per atom type:
    std::vector<int> natoms(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        natoms[i] = std::count(data.atnum.begin(), data.atnum.end(), atoms[i]);
    }
    // Create Dalton mol file:
    using namespace Chem::Periodic_table;

    Stdutils::Format<double> fix;
    fix.fixed().width(12).precision(6);

    std::cout << "Angstrom Atomtypes=" << atoms.size() << '\n';
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        std::cout << "Charge=" << atoms[i] << " Atoms=" << natoms[i] << '\n';
        int iter = 1;
        for (std::size_t j = 0; j < data.atnum.size(); ++j) {
            if (data.atnum[j] == atoms[i]) {
                std::cout << get_atomic_symbol(atoms[i]) << iter << '\t'
                          << fix(data.xyz(j, 0)) << ' ' << fix(data.xyz(j, 1))
                          << ' ' << fix(data.xyz(j, 2)) << '\n';
                ++iter;
            }
        }
    }
}
