// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_IO_H
#define CHEM_IO_H

#include <chem/element.h>
#include <chem/traits.h>
#include <numlib/matrix.h>
#include <iostream>
#include <vector>

namespace Chem {

// Read chemical XYZ file format.
void read_xyz_format(std::istream& from,
                     std::vector<Element>& atoms,
                     Numlib::Mat<double>& xyz,
                     std::string& title);

// Read chemical Z matrix format.
void read_zmat_format(std::istream& from,
                      std::vector<Element>& atoms,
                      Numlib::Vec<double>& distances,
                      Numlib::Vec<double>& angles,
                      Numlib::Vec<double>& dihedrals,
                      Numlib::Vec<Index>& bond_connect,
                      Numlib::Vec<Index>& angle_connect,
                      Numlib::Vec<Index>& dihedral_connect);

// Read molecular formula.
void read_mol_formula(std::istream& from, std::vector<Mol_formula>& formula);

// Print chemical XYZ file format.
void print_xyz_format(std::ostream& to,
                      const std::vector<Element>& atoms,
                      const Numlib::Mat<double>& xyz,
                      const std::string& title);

// Print chemical Z matrix format.
void print_zmat_format(std::ostream& to,
                       const std::vector<Element>& atoms,
                       const Numlib::Vec<double>& distances,
                       const Numlib::Vec<double>& angles,
                       const Numlib::Vec<double>& dihedrals,
                       const Numlib::Vec<Index>& bond_connect,
                       const Numlib::Vec<Index>& angle_connect,
                       const Numlib::Vec<Index>& dihedral_connect);

// Print molecular electronic states.
void print_spin_orbit_states(std::ostream& to,
                             const Numlib::Vec<int>& so_degen,
                             const Numlib::Vec<double>& so_energy);

// Print molecular geometry.
void print_geometry(std::ostream& to,
                    const std::vector<Element>& atoms,
                    const Numlib::Mat<double>& xyz,
                    const std::string& unit = "angstrom");

// Print atomic masses.
void print_atomic_masses(std::ostream& to, const std::vector<Element>& atoms);

// Print center-of-mass coordinates.
void print_center_of_mass(std::ostream& to, const Numlib::Vec<double>& com);

// Print principal moments.
void print_principal_moments(std::ostream& to,
                             const Numlib::Vec<double>& pmom,
                             const Numlib::Mat<double>& paxis);

// Print rotational constants.
void print_rot_constants(std::ostream& to,
                         int sigma,
                         const std::string& symm,
                         const Numlib::Vec<double>& rotc);

} // namespace Chem

#endif // CHEM_IO_H

