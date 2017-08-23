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

#ifndef CHEM_MOLECULE_H
#define CHEM_MOLECULE_H

#include <chem/element.h>
#include <chem/molrot.h>
#include <chem/molvib.h>
#include <chem/torsion.h>
#include <chem/zmatrix.h>
#include <armadillo>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Error reporting:

struct Mol_error : std::runtime_error {
    Mol_error(std::string s) : std::runtime_error(s) {}
};

//
// Class for holding molecule data.
//
class Molecule {
public:
    Molecule()
    {
        zmat = std::make_unique<Zmatrix>(atoms, xyz);
        rot  = std::make_unique<Molrot>(atoms, xyz);
        vib  = std::make_unique<Molvib>();
        tor  = std::make_unique<Torsion>(*rot);
    }

    Molecule(std::istream& from,
             std::ostream& to       = std::cout,
             const std::string& key = "Molecule",
             bool verbose           = false)
    {
        init(from, to, key, verbose);
    }

    Molecule(const Molecule& mol);

    ~Molecule() {}

    // Calculate total molecular mass.
    double tot_mass() const { return rot->tot_mass(); }

    const std::string get_title() const { return title; }
    const std::vector<Element>& get_atoms() const { return atoms; }
    const arma::mat& get_xyz() const { return xyz; }
    const arma::vec& get_elec_state() const { return elec_state; }

    Zmatrix& get_zmat() const { return *zmat; }
    Molrot& get_rot() const { return *rot; }
    Molvib& get_vib() const { return *vib; }
    Torsion& get_tor() const { return *tor; }

    int get_charge() const { return charge; }
    double get_elec_energy() const { return elec_energy; }

    void set_xyz(const arma::mat& xyz_);
    void set_charge(const int charge_) { charge = charge_; }
    void set_elec_energy(const double energy) { elec_energy = energy; }

    void print_data(std::ostream& to, const std::string& key) const;

private:
    void init(std::istream& from,
              std::ostream& to,
              const std::string& key,
              bool verbose);

    std::string title;      // molecule information
    std::string geom_unit;  // units for geometry

    std::vector<Element> atoms;  // atoms in molecule
    arma::mat xyz;               // cartesian coordinates
    arma::vec elec_state;        // electronic state

    int charge;          // molecular charge
    double elec_energy;  // electronic ground-state energy [Hartree]

    std::unique_ptr<Zmatrix> zmat;
    std::unique_ptr<Molrot> rot;
    std::unique_ptr<Molvib> vib;
    std::unique_ptr<Torsion> tor;
};

inline void Molecule::set_xyz(const arma::mat& xyz_)
{
    // Note: Moment of inertia for torsional modes is currently not updated
    // when the geometry is changed. This may change in future versions.
    //
    // TODO (stigrs@gmail.com) Check if torsional modes need to be updated.
    xyz = xyz_;
    zmat->build_zmat();
}

#endif  // CHEM_MOLECULE_H
