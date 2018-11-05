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
#include <chem/impl/io_support.h>
#include <numlib/matrix.h>
#include <vector>

namespace Chem {

// Class for handling Z matrices.
//
class Zmatrix {
public:
    Zmatrix() = delete;

    // Construct Z matrix from atoms and Cartesian coordinates.
    Zmatrix(std::vector<Element>& atoms_, Numlib::Mat<double>& xyz_);

    // Copy semantics:
    Zmatrix(const Zmatrix&) = default;
    Zmatrix& operator=(const Zmatrix&) = default;

    // Move semantics:
    Zmatrix(Zmatrix&&) = default;
    Zmatrix& operator=(Zmatrix&&) = default;

    ~Zmatrix() = default;

    // Get bond distance.
    double get_distance(int index) const;

    // Get bond angle.
    double get_angle(int index) const;

    // Get dihedral angle.
    double get_dihedral(int index) const;

    // Get connectivities.
    std::vector<Numlib::Vec<int>> get_connectivities() const;

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
    void print(std::ostream& to = std::cout) const;

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
    int find_nearest_atom(const Numlib::Vec<double>& dist) const;

    // Find new connection.
    int find_new_connection(const Numlib::Vec<int>& iatms,
                            const Numlib::Vec<int>& connectivity) const;

    // Calculate position of another atom based on internal coordinates. The
    // code is based on the qcl code written by Ben Albrecht released under
    // the MIT license.
    Numlib::Vec<double> calc_position(int i) const;

private:
    std::vector<Element>& atoms;
    Numlib::Mat<double>& xyz;

    Numlib::Vec<double> distances;
    Numlib::Vec<double> angles;
    Numlib::Vec<double> dihedrals;

    Numlib::Vec<int> bond_connect;
    Numlib::Vec<int> angle_connect;
    Numlib::Vec<int> dihedral_connect;
};

inline Zmatrix::Zmatrix(std::vector<Element>& atoms_, Numlib::Mat<double>& xyz_)
    : atoms(atoms_), xyz(xyz_)
{
    if (!atoms.empty()) {
        distances.resize(atoms.size());
        angles.resize(atoms.size());
        dihedrals.resize(atoms.size());

        bond_connect.resize(atoms.size());
        angle_connect.resize(atoms.size());
        dihedral_connect.resize(atoms.size());

        build_zmat();
    }
}

inline double Zmatrix::get_distance(int index) const
{
    double res = 0.0;
    if (atoms.size() > 1) {
        res = distances(index);
    }
    return res;
}

inline double Zmatrix::get_angle(int index) const
{
    double res = 0.0;
    if (atoms.size() > 2) {
        res = angles(index);
    }
    return res;
}

inline double Zmatrix::get_dihedral(int index) const
{
    double res = 0.0;
    if (atoms.size() > 3) {
        res = dihedrals(index);
    }
    return res;
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

inline void Zmatrix::load(std::istream& from)
{
    Impl::read_zmat_format(from, atoms, distances, angles, dihedrals,
                           bond_connect, angle_connect, dihedral_connect);
    build_xyz();
}

inline void Zmatrix::print(std::ostream& to) const
{
    Impl::print_zmat_format(to, atoms, distances, angles, dihedrals,
                            bond_connect, angle_connect, dihedral_connect);
}

} // namespace Chem

#endif // CHEM_ZMATRIX_H

