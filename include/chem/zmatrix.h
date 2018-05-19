////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CHEM_ZMATRIX_H
#define CHEM_ZMATRIX_H

#include <chem/element.h>
#include <srs/array.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Error reporting:

struct Zmatrix_error : std::runtime_error {
    Zmatrix_error(std::string s) : std::runtime_error(s) {}
};

// Forward declarations to allow friend declarations:

class Molecule;

//
// Class for handling Z matrices.
//
class Zmatrix {
public:
    Zmatrix(std::vector<Element>& atoms_, srs::dmatrix& xyz_);

    Zmatrix(const Zmatrix& zmat);

    // Get bond distance.
    double get_distance(int index) const;

    // Get bond angle.
    double get_angle(int index) const;

    // Get dihedral angle.
    double get_dihedral(int index) const;

    // Get connectivities.
    std::vector<srs::ivector> get_connectivities() const;

    // Set bond distance.
    void set_distance(int index, double value);

    // Set bond angle.
    void set_angle(int index, double value);

    // Set dihedral angle.
    void set_dihedral(int index, double value);

    // Rotate moiety around a given torsional bond.
    void rotate_moiety(const std::vector<int>& moiety, double value);

    // Load molecular coordinates in Z matrix format.
    void load(std::istream& from);

    // Print Z matrix.
    void print(std::ostream& to);

    friend class Molecule;

protected:
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
    int find_nearest_atom(const srs::dvector& dist) const;

    // Find new connection.
    int find_new_connection(const srs::ivector& iatms,
                            const srs::ivector& connectivity) const;

    // Calculate position of another atom based on internal coordinates. The
    // code is based on the qcl code written by Ben Albrecht released under
    // the MIT license.
    srs::dvector calc_position(int i) const;

private:
    std::vector<Element>& atoms;
    srs::dmatrix& xyz;

    srs::dvector distances;
    srs::dvector angles;
    srs::dvector dihedrals;

    srs::ivector bond_connect;
    srs::ivector angle_connect;
    srs::ivector dihedral_connect;
};

inline Zmatrix::Zmatrix(std::vector<Element>& atoms_, srs::dmatrix& xyz_)
    : atoms(atoms_), xyz(xyz_)
{
    distances.resize(atoms.size(), 0.0);
    angles.resize(atoms.size(), 0.0);
    dihedrals.resize(atoms.size(), 0.0);

    bond_connect.resize(atoms.size(), 0);
    angle_connect.resize(atoms.size(), 0);
    dihedral_connect.resize(atoms.size(), 0);

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

#endif  // CHEM_ZMATRIX_H
