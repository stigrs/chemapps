// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_MOLECULE_IO_H
#define CHEM_MOLECULE_IO_H

#include <chem/element.h>
#include <numlib/matrix.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

namespace Impl {

    // Read chemical XYZ file format.
    void read_xyz_format(std::istream& from,
                         std::vector<Chem::Element>& atoms,
                         Numlib::Mat<double>& xyz,
                         std::string& title);

    // Read chemical Z matrix format.
    void read_zmat_format(std::istream& from,
                          std::vector<Chem::Element>& atoms,
                          Numlib::Vec<double>& distances,
                          Numlib::Vec<double>& angles,
                          Numlib::Vec<double>& dihedrals,
                          Numlib::Vec<int>& bond_connect,
                          Numlib::Vec<int>& angle_connect,
                          Numlib::Vec<int>& dihedral_connect);

    // Print molecular electronic states.
    void print_elec_states(const Numlib::Vec<double>& elec_state);

    // Print chemical Z matrix format.
    void print_zmat_format(const std::vector<Chem::Element>& atoms,
                           const Numlib::Vec<double>& distances,
                           const Numlib::Vec<double>& angles,
                           const Numlib::Vec<double>& dihedrals,
                           const Numlib::Vec<int>& bond_connect,
                           const Numlib::Vec<int>& angle_connect,
                           const Numlib::Vec<int>& dihedral_connect);

}  // namespace Impl

}  // namespace Chem

#endif  // CHEM_MOLECULE_IO_H
