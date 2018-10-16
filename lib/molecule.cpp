// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/molecule.h>

Chem::Molecule::Molecule(std::istream& from,
                         const std::string& key,
                         bool verbose)
    : elec(from, key),
      geom(from, key),
      rot(from, key, geom),
      vib(from, key, geom, rot)
{
    if (verbose) {
        std::cout << "verbose\n";
    }
}

