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

#ifndef CHEM_MOLECULE_IO_H
#define CHEM_MOLECULE_IO_H

#include <chem/element.h>
#include <chem/traits.h>
#include <numlib/matrix.h>
#include <iostream>
#include <vector>

namespace Chem {

namespace Impl {

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
                          Numlib::Vec<int>& bond_connect,
                          Numlib::Vec<int>& angle_connect,
                          Numlib::Vec<int>& dihedral_connect);

    // Read molecular formula.
    void read_mol_formula(std::istream& from,
                          std::vector<Mol_formula>& formula);

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
                           const Numlib::Vec<int>& bond_connect,
                           const Numlib::Vec<int>& angle_connect,
                           const Numlib::Vec<int>& dihedral_connect);

    // Print molecular electronic states.
    void print_elec_states(std::ostream& to,
                           const Numlib::Vec<double>& elec_state);

    // Print molecular geometry.
    void print_geometry(std::ostream& to,
                        const std::vector<Element>& atoms,
                        const Numlib::Mat<double>& xyz,
                        const std::string& unit = "angstrom");

    // Print atomic masses.
    void print_atomic_masses(std::ostream& to,
                             const std::vector<Element>& atoms);

} // namespace Impl

} // namespace Chem

#endif  // CHEM_MOLECULE_IO_H

