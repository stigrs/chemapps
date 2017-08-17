/**
   @file zmatrix.h

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

#ifndef CHEM_ZMATRIX_H
#define CHEM_ZMATRIX_H

#include <chem/element.h>
#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------

// Error reporting:

struct Zmatrix_error : std::runtime_error {
    Zmatrix_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

// Forward declarations to allow friend declarations:

class Molecule;

//-----------------------------------------------------------------------------

class Zmatrix {
public:
    Zmatrix(std::vector<Element>& atoms_, arma::mat& xyz_);

    ~Zmatrix() {}

    /// Get bond distance.
    double get_distance(int index) const;

    /// Get bond angle.
    double get_angle(int index) const;

    /// Get dihedral angle.
    double get_dihedral(int index) const;

    /// Get connectivities.
    std::vector<arma::ivec> get_connectivities() const;

    /// Set bond distance.
    void set_distance(int index, double value);

    /// Set bond angle.
    void set_angle(int index, double value);

    /// Set dihedral angle.
    void set_dihedral(int index, double value);

    /// Rotate moiety around a given torsional bond.
    void rotate_moiety(const arma::ivec& moiety, double value);

    /// Load molecular coordinates in Z matrix format.
    void load(std::istream& from);

    /// Print Z matrix.
    void print(std::ostream& to);

    friend class Molecule;

protected:
    /**
       Build Z matrix from Cartesian coordinates. The code is based on the
       qcl code written by Ben Albrecht released under the MIT license.

       @note It is assumed that bonded atoms are closer than non-bonded atoms.
       This may not work well for transition states and molecular complexes.
    */
    void build_zmat();

    /**
       Convert Z matrix to Cartesian coordinates. The code is based on the
       qcl code written by Ben Albrecht released under the MIT license.
    */
    void build_xyz();

    /// Find index to nearest atom.
    int find_nearest_atom(const arma::rowvec& dist) const;

    /// Find new connection.
    int find_new_connection(const arma::ivec& iatms,
                            const arma::ivec& connectivity) const;

    /**
       Calculate position of another atom based on internal coordinates. The
       code is based on the qcl code written by Ben Albrecht released under
       the MIT license.
    */
    arma::rowvec calc_position(int i) const;

private:
    std::vector<Element>& atoms;
    arma::mat& xyz;

    arma::vec distances;
    arma::vec angles;
    arma::vec dihedrals;

    arma::ivec bond_connect;
    arma::ivec angle_connect;
    arma::ivec dihedral_connect;
};

inline Zmatrix::Zmatrix(std::vector<Element>& atoms_, arma::mat& xyz_)
    : atoms(atoms_), xyz(xyz_)
{
    distances = arma::zeros<arma::vec>(atoms.size());
    angles    = arma::zeros<arma::vec>(atoms.size());
    dihedrals = arma::zeros<arma::vec>(atoms.size());

    bond_connect     = arma::zeros<arma::ivec>(atoms.size());
    angle_connect    = arma::zeros<arma::ivec>(atoms.size());
    dihedral_connect = arma::zeros<arma::ivec>(atoms.size());

    build_zmat();
}

inline double Zmatrix::get_distance(int index) const
{
    if (atoms.size() > 1) {
        return distances(index);
    }
    else {
        return 0.0;
    }
}

inline double Zmatrix::get_angle(int index) const
{
    if (atoms.size() > 2) {
        return angles(index);
    }
    else {
        return 0.0;
    }
}

inline double Zmatrix::get_dihedral(int index) const
{
    if (atoms.size() > 3) {
        return dihedrals(index);
    }
    else {
        return 0.0;
    }
}

inline void Zmatrix::set_distance(int index, double value)
{
    if (atoms.size() > 1) {
        distances(index) = value;
        build_xyz();
    }
}

inline void Zmatrix::set_angle(int index, double value)
{
    if (atoms.size() > 2) {
        angles(index) = value;
        build_xyz();
    }
}

inline void Zmatrix::set_dihedral(int index, double value)
{
    if (atoms.size() > 3) {
        dihedrals(index) = value;
        build_xyz();
    }
}

#endif /* CHEM_ZMATRIX_H */
