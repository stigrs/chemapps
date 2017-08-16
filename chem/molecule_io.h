/**
   @file molecule_io.h

   This file is part of ChemApps - A C++ Chemistry Toolkit

   Copyright (C) 2016-2017  Stig Rune Sellevag

   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CHEM_MOLECULE_IO_H
#define CHEM_MOLECULE_IO_H

#include <chem/element.h>
#include <armadillo>
#include <iostream>
#include <vector>

namespace chem {

    /// Read chemical XYZ file format.
    void read_xyz_format(std::istream& from,
                         std::vector<Element>& atoms,
                         arma::mat& xyz,
                         std::string& title);

    /// Read chemical Z matrix format.
    void read_zmat_format(std::istream& from,
                          std::vector<Element>& atoms,
                          arma::vec& distances,
                          arma::vec& angles,
                          arma::vec& dihedrals,
                          arma::ivec& bond_connect,
                          arma::ivec& angle_connect,
                          arma::ivec& dihedral_connect);

    /// Print chemical XYZ file format.
    void print_xyz_format(std::ostream& to,
                          const std::vector<Element>& atoms,
                          const arma::mat& xyz,
                          const std::string& title);

    /// Print chemical Z matrix format.
    void print_zmat_format(std::ostream& to,
                           std::vector<Element>& atoms,
                           arma::vec& distances,
                           arma::vec& angles,
                           arma::vec& dihedrals,
                           arma::ivec& bond_connect,
                           arma::ivec& angle_connect,
                           arma::ivec& dihedral_connect);

    /// Print molecular electronic states.
    void print_elec_states(std::ostream& to, const arma::vec& elec_state);

    /// Print molecular geometry.
    void print_geometry(std::ostream& to,
                        const std::vector<Element>& atoms,
                        const arma::mat& xyz,
                        const std::string& unit = "angstrom");

    /// Print atomic masses.
    void print_atomic_masses(std::ostream& to,
                             const std::vector<Element>& atoms);
}  // chem::

#endif /* CHEM_MOLECULE_IO_H */
