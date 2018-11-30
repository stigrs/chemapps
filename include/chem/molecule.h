// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_MOLECULE_H
#define CHEM_MOLECULE_H

#include <chem/electronic.h>
#include <chem/geometry.h>
#include <chem/rotation.h>
#include <chem/vibration.h>
#include <chem/torsion.h>
#include <chem/traits.h>
#include <iostream>
#include <string>

namespace Chem {

// Class for holding molecule objects.
//
class Molecule {
public:
    Molecule() = default;

    Molecule(std::istream& from,
             std::ostream& to = std::cout,
             const std::string& key = "Molecule",
             bool verbose = false);

    // Copy semantics:
    Molecule(const Molecule&) = default;
    Molecule& operator=(const Molecule&) = default;

    // Move semantics:
    Molecule(Molecule&&) = default;
    Molecule& operator=(Molecule&&) = default;

    ~Molecule() = default;

    // Get information string for molecule.
    std::string title() const { return geom_.title(); }

    // Get number of atoms.
    auto num_atoms() const { return geom_.atoms().size(); }

    // Get atoms in molecule.
    const auto& atoms() const { return geom_.atoms(); }

    // Get total molecular mass.
    double tot_mass() const;

    // Get electronic properties.
    const auto& elec() const { return elec_; }

    // Modify electronic properties.
    auto& elec() { return elec_; }

    // Get molecular structure.
    Mol_type structure() const;

    // Get Cartesian coordinates.
    const auto& get_xyz() const { return geom_.get_xyz(); }

    // Set molecular geometry.
    void set_xyz(const Numlib::Mat<double>& x);

    // Get molecular geometry.
    const auto& geom() const { return geom_; }

    // Modify molecular geometry.
    auto& geom() { return geom_; }

    // Get molecular rotation object.
    const auto& rot() const { return rot_; }

    // Modify molecular rotations.
    auto& rot() { return rot_; }

    // Get molecular vibrations object.
    const auto& vib() const { return vib_; }

    // Modify molecular vibrations.
    auto& vib() { return vib_; }

    // Get internal torsions object.
    const auto& tor() const { return tor_; }

    // Modify internal torsions.
    auto& tor() { return tor_; }

private:
    Electronic elec_; // electronic properties
    Geometry geom_;   // molecular geometry
    Rotation rot_;    // molecular rotations
    Vibration vib_;   // molecular vibrations
    Torsion tor_;     // internal torsional modes
};

inline Mol_type Molecule::structure() const
{
    Mol_type res;
    if (num_atoms() == 1) {
        res = atom;
    }
    else if (rot_.constants().size() == 1) {
        res = linear;
    }
    else {
        res = nonlinear;
    }
    return res;
}

inline void Molecule::set_xyz(const Numlib::Mat<double>& x)
{
    Assert::dynamic(Numlib::same_extents(geom_.get_xyz(), x),
                    "bad size of Cartesian coordinates");
    geom_.set_xyz(x);
    rot_.set(x);
    vib_ = Vibration(); // geometry has changed; frequencies no longer valid
    tor_.set(x, rot_.principal_axes(), rot_.principal_moments());
}

} // namespace Chem

#endif // CHEM_MOLECULE_H

