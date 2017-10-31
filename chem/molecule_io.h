//////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef CHEM_MOLECULE_IO_H
#define CHEM_MOLECULE_IO_H

#include <chem/element.h>
#include <srs/array.h>
#include <iostream>
#include <vector>

namespace chem {

// Read chemical XYZ file format.
void read_xyz_format(std::istream& from,
                     std::vector<Element>& atoms,
                     srs::dmatrix& xyz,
                     std::string& title);

// Read chemical Z matrix format.
void read_zmat_format(std::istream& from,
                      std::vector<Element>& atoms,
                      srs::dvector& distances,
                      srs::dvector& angles,
                      srs::dvector& dihedrals,
                      srs::ivector& bond_connect,
                      srs::ivector& angle_connect,
                      srs::ivector& dihedral_connect);

// Print chemical XYZ file format.
void print_xyz_format(std::ostream& to,
                      const std::vector<Element>& atoms,
                      const srs::dmatrix& xyz,
                      const std::string& title);

// Print chemical Z matrix format.
void print_zmat_format(std::ostream& to,
                       std::vector<Element>& atoms,
                       srs::dvector& distances,
                       srs::dvector& angles,
                       srs::dvector& dihedrals,
                       srs::ivector& bond_connect,
                       srs::ivector& angle_connect,
                       srs::ivector& dihedral_connect);

// Print molecular electronic states.
void print_elec_states(std::ostream& to, const srs::dvector& elec_state);

// Print molecular geometry.
void print_geometry(std::ostream& to,
                    const std::vector<Element>& atoms,
                    const srs::dmatrix& xyz,
                    const std::string& unit = "angstrom");

// Print atomic masses.
void print_atomic_masses(std::ostream& to, const std::vector<Element>& atoms);

}  // namespace chem

#endif  // CHEM_MOLECULE_IO_H
