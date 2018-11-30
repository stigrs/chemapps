// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_WHIRAB_H
#define CHEM_WHIRAB_H

#include <chem/molecule.h>

namespace Chem {

// Provides Whitten-Rabinovitch approximations.
//
namespace Whirab {

    // Whitten-Rabinovitch correction.
    double a_corr(const Molecule& mol, double e_barrier);

    // Vibrational density of states.
    double vibr_density_states(const Molecule& mol, double e_barrier);

} // namespace Whirab

} // namespace Chem

#endif // CHEM_WHIRAB_H

