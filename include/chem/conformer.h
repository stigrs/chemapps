// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_CONFORMER_H
#define CHEM_CONFORMER_H

#include <chem/element.h>
#include <numlib/matrix.h>
#include <vector>

namespace Chem {

// Simple structure for storing conformers.
//
struct Conformer {
    Conformer() = default;

    explicit Conformer(double e, const Numlib::Mat<double>& x)
        : energy{e}, atoms{0}, xyz(x), iter{0}
    {
    }

    Conformer(double e,
              const std::vector<Element>& at,
              const Numlib::Mat<double>& x)
        : energy{e}, atoms(at), xyz(x), iter{0}
    {
    }

    // Copy semantics:
    Conformer(const Conformer&) = default;
    Conformer& operator=(const Conformer&) = default;

    // Move semantics:
    Conformer(Conformer&&) = default;
    Conformer& operator=(Conformer&&) = default;

    ~Conformer() = default;

    // Compare conformers by energy.
    bool operator<(const Conformer& c) const { return energy < c.energy; }

    // Compare conformers by energy.
    bool operator>(const Conformer& c) const { return energy > c.energy; }

    double energy;              // conformer energy
    std::vector<Element> atoms; // atoms
    Numlib::Mat<double> xyz;    // Cartesian coordinates
    int iter;                   // iterator to be used for counting
};

} // namespace Chem

#endif // CHEM_CONFORMER_H

