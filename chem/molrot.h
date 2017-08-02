/**
   @file molrot.h
   
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

#ifndef CHEM_MOLROT_H
#define CHEM_MOLROT_H

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <armadillo>
#include <chem/element.h>
#include <chem/math.h>

//-----------------------------------------------------------------------------

// Error reporting:

struct Molrot_error : std::runtime_error {
    Molrot_error(std::string s) : std::runtime_error(s) { }
};

//-----------------------------------------------------------------------------

// Forward declarations to allow friend declarations:

class Imom_tor;

//-----------------------------------------------------------------------------

/// Class for handling molecular rotations.
class Molrot {
public:
    Molrot(std::vector<Element>& atoms_, arma::mat& xyz_)
        : atoms(atoms_), xyz(xyz_)
    { }

    Molrot(std::istream& from, 
           const std::string& key,
           std::vector<Element>& atoms_,
           arma::mat& xyz_);

    ~Molrot() { }

    /// Perform rotational analysis.
    void analysis(std::ostream& to = std::cout);

    /// Compute rotational constants.
    arma::vec3 constants();

    /// Return rotational symmetry number.
    double get_sigma() const { return sigma; }

    /// Return rotational symmetry.
    std::string symmetry();

    friend class Imom_tor;
    
protected:
    /// Initialize.
    void init(std::istream& from, const std::string& key);

    /// Move geometry to center of mass.
    void move_to_com();

    /// Compute center of mass coordinates.
    arma::vec3 center_of_mass() const;

    /// Compute principal moments of inertia.
    void principal_moments();

    /**
       Rotate to principal axes.

       @note This coordinate system is not the same as Gaussian's 
       standard orientation.
    */
    void rotate_to_principal_axes();

    /// Print center of mass.
    void print_center_of_mass(std::ostream& to = std::cout) const;

    /// Print principal moments and axes.
    void print_principal_moments(std::ostream& to = std::cout) const;

    /// Print rotational constants.
    void print_constants(std::ostream& to = std::cout);

    std::vector<Element>& atoms;
    arma::mat& xyz;

    arma::vec3  pmom;
    arma::mat33 paxis;

    double sigma;          
    bool   aligned = false;
 
    Molrot& operator=(const Molrot& rot); // no assignments
};

inline Molrot::Molrot(std::istream& from,
                      const std::string& key,
                      std::vector<Element>& atoms_,
                      arma::mat& xyz_)
    : atoms(atoms_), xyz(xyz_)
{
    init(from, key);
}

inline void Molrot::move_to_com()
{
    arma::vec3 com = center_of_mass();
    chem::translate(xyz, -com(0), -com(1), -com(2));
}

#endif /* CHEM_MOLROT_H */

