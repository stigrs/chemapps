// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEMLIB_TORSION_H
#define CHEMLIB_TORSION_H

#include <chem/impl/geometry.h>
#include <chem/impl/rotation.h>

namespace Chem {

namespace Impl {

    // Class for handling torsional modes in molecules using the CT-Cw scheme.
    //
    // Algorithm:
    // ----------
    // The reduced moment of inertia of a symmetrical or unsymmetrical rotating
    // top attached to a rigid frame is calculated according to eq 1 in
    // Pitzer, K. S. J. Chem. Phys. 1946, vol. 14, pp. 239-243, also known as
    // the curvilinear (C) scheme.
    //
    // The CT-Cw scheme is described in the following paper: Chuang, Y.-Y.;
    //  Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.
    //
    // Note:
    // -----
    // Atoms specifying the rotational axis must not be included in the list
    // of atoms specifying the rotating top.
    //
    class Torsion {
    public:
        Torsion() = delete;

        Torsion
    };

} // Molecule_impl

} // Chemlib

#endif // CHEMLIB_TORSION_H
