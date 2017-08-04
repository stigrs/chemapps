/**
   @file molecule.h
   
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

#ifndef CHEM_MOLECULE_H
#define CHEM_MOLECULE_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <memory>
#include <armadillo>
#include <chem/element.h>
#include <chem/zmatrix.h>
#include <chem/molrot.h>
#include <chem/molvib.h>
#include <chem/torsion.h>

//----------------------------------------------------------------------------

// Error reporting:

struct Mol_error : std::runtime_error {
    Mol_error(std::string s) : std::runtime_error(s) { }
};

//----------------------------------------------------------------------------

/// Class for holding molecule data.
class Molecule {
public:
    Molecule() 
    { 
        zmat = std::make_shared<Zmatrix>(atoms, xyz);
        rot  = std::make_shared<Molrot>(atoms, xyz);
        vib  = std::make_shared<Molvib>();
        tor  = std::make_shared<Torsion>(*rot);
    }

    Molecule(std::istream& from, 
             std::ostream& to, 
             const std::string& key,
             bool verbose = false)
    {
        init(from, to, key, verbose);
    }

    ~Molecule() { }

    std::string get_title() const { return title; }

    const std::vector<Element>& get_atoms()      const { return atoms; }
    const arma::mat&            get_xyz()        const { return xyz; }
    const arma::vec&            get_elec_state() const { return elec_state; }

    std::shared_ptr<Zmatrix> get_zmat() const { return zmat; }
    std::shared_ptr<Molrot>  get_rot()  const { return rot; }
    std::shared_ptr<Molvib>  get_vib()  const { return vib; }
    std::shared_ptr<Torsion> get_tor()  const { return tor; }

    int    get_charge()      const { return charge; }
    double get_elec_energy() const { return elec_energy; }

    void set_xyz(const arma::mat& xyz_)       { xyz = xyz_; }
    void set_charge(const int charge_)        { charge = charge_; }
    void set_elec_energy(const double energy) { elec_energy = energy; }

    void print_data(std::ostream& to, const std::string& key) const;

private:
    void init(std::istream& from, 
              std::ostream& to, 
              const std::string& key,
              bool verbose);
    
    std::string title;     ///< molecule information
    std::string geom_unit; ///< units for geometry

    std::vector<Element> atoms;      ///< atoms in molecule 
    arma::mat            xyz;        ///< cartesian coordinates
    arma::vec            elec_state; ///< electronic state

    int    charge;      ///< molecular charge
    double elec_energy; ///< electronic ground-state energy [Hartree]
   
    std::shared_ptr<Zmatrix> zmat;
    std::shared_ptr<Molrot>  rot;
    std::shared_ptr<Molvib>  vib;
    std::shared_ptr<Torsion> tor;
};

#endif /* CHEM_MOLECULE_H */

