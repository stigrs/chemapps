// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_GEOMETRY_H
#define CHEM_GEOMETRY_H

#include <chem/element.h>
#include <chem/io.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

// Class for handling molecular geometries.
//
class Geometry {
public:
    Geometry() = default;

    Geometry(std::istream& from, const std::string& key);

    // Copy semantics:
    Geometry(const Geometry&) = default;
    Geometry& operator=(const Geometry&) = default;

    // Move semantics:
    Geometry(Geometry&&) = default;
    Geometry& operator=(Geometry&&) = default;

    ~Geometry() = default;

    // Get title string:
    std::string title() const { return info; }

    // Get atoms.
    const auto& atoms() const { return atms; }

    // Get Cartesian coordinates.
    const auto& get_xyz() const { return xyz; }

    // Get bond distance.
    double get_distance(Index index) const;

    // Get bond angle.
    double get_angle(Index index) const;

    // Get dihedral angle.
    double get_dihedral(Index index) const;

    // Get connectivities.
    std::vector<Numlib::Vec<Index>> get_connectivities() const;

    // Set Cartesian coordinates.
    void set_xyz(const Numlib::Mat<double>& x);

    // Set bond distance.
    void set_distance(Index index, double value);

    // Set bond angle.
    void set_angle(Index index, double value);

    // Set dihedral angle.
    void set_dihedral(Index index, double value);

    // Rotate moiety around a given torsional bond.
    void rotate_moiety(const std::vector<Index>& moiety, double value);

    // Load molecular coordinates in XYZ format.
    void load_xyz(std::istream& from);

    // Load molecular coordinates in Z matrix format.
    void load_zmat(std::istream& from);

    // Print molecular coordinates in XYZ format.
    void print_xyz(std::ostream& to = std::cout) const;

    // Print molecular coordinates in Z matrix format.
    void print_zmat(std::ostream& to = std::cout) const;

private:
    // Build Z matrix from Cartesian coordinates. The code is based on the
    // qcl code written by Ben Albrecht released under the MIT license.
    //
    // Note: It is assumed that bonded atoms are closer than non-bonded atoms.
    // This may not work well for transition states and molecular complexes.
    void build_zmat();

    // Convert Z matrix to Cartesian coordinates. The code is based on the
    // qcl code written by Ben Albrecht released under the MIT license.
    void build_xyz();

    // Find index to nearest atom.
    Index find_nearest_atom(const Numlib::Vec<double>& dist) const;

    // Find new connection.
    Index find_new_connection(const Numlib::Vec<Index>& iatms,
                              const Numlib::Vec<Index>& connectivity) const;

    // Calculate position of another atom based on internal coordinates. The
    // code is based on the qcl code written by Ben Albrecht released under
    // the MIT license.
    Numlib::Vec<double> calc_position(Index i) const;

    std::vector<Element> atms;
    Numlib::Mat<double> xyz;

    Numlib::Vec<double> distances;
    Numlib::Vec<double> angles;
    Numlib::Vec<double> dihedrals;

    Numlib::Vec<Index> bond_connect;
    Numlib::Vec<Index> angle_connect;
    Numlib::Vec<Index> dihedral_connect;

    std::string info;
};

inline double Geometry::get_distance(Index index) const
{
    double res = 0.0;
    if (atms.size() > 1) {
        res = distances(index);
    }
    return res;
}

inline double Geometry::get_angle(Index index) const
{
    double res = 0.0;
    if (atms.size() > 2) {
        res = angles(index);
    }
    return res;
}

inline double Geometry::get_dihedral(Index index) const
{
    double res = 0.0;
    if (atms.size() > 3) {
        res = dihedrals(index);
    }
    return res;
}

inline void Geometry::set_xyz(const Numlib::Mat<double>& x)
{
    Assert::dynamic(Numlib::same_extents(xyz, x),
                    "bad size of Cartesian coordinates");
    xyz = x;
    build_zmat();
}

inline void Geometry::set_distance(Index index, double value)
{
    if (atms.size() > 1) {
        distances(index) = value;
        build_xyz();
    }
}

inline void Geometry::set_angle(Index index, double value)
{
    if (atms.size() > 2) {
        angles(index) = value;
        build_xyz();
    }
}

inline void Geometry::set_dihedral(Index index, double value)
{
    if (atms.size() > 3) {
        dihedrals(index) = value;
        build_xyz();
    }
}

inline void Geometry::load_xyz(std::istream& from)
{
    read_xyz_format(from, atms, xyz, info);
    build_zmat();
}

inline void Geometry::load_zmat(std::istream& from)
{
    read_zmat_format(from, atms, distances, angles, dihedrals, bond_connect,
                     angle_connect, dihedral_connect);
    build_xyz();
}

inline void Geometry::print_xyz(std::ostream& to) const
{
    print_xyz_format(to, atms, xyz, info);
}

inline void Geometry::print_zmat(std::ostream& to) const
{
    print_zmat_format(to, atms, distances, angles, dihedrals, bond_connect,
                      angle_connect, dihedral_connect);
}

} // namespace Chem

#endif // CHEM_GEOMETRY_H
