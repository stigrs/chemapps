/**
   @file torsion.cpp
   
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

#include <map>
#include <cmath>
#include <chem/torsion.h>
#include <chem/input.h>
#include <chem/utils.h>


void Torsion::init(std::istream& from, const std::string& key)
{
    typedef std::map<std::string, Input>::iterator       Input_iter;
    typedef std::map<std::string, Input>::const_iterator Cinput_iter;

    // Read input data:

    arma::uvec rot_axis_def  = arma::zeros<arma::uvec>(2);
    arma::uvec rot_top_def   = arma::zeros<arma::uvec>(0);
    arma::uvec sigma_tor_def = arma::zeros<arma::uvec>(0);
    arma::vec  rmi_tor_def   = arma::zeros<arma::vec>(0);
    arma::vec  pot_tor_def   = arma::zeros<arma::vec>(0);
    arma::vec  freq_tor_def  = arma::zeros<arma::vec>(0);
    
    std::map<std::string, Input> input_data;
    input_data["rot_axis"]  = Input(rot_axis, rot_axis_def);
    input_data["rot_top"]   = Input(rot_top, rot_top_def);
    input_data["sigma_tor"] = Input(sigma_tor, sigma_tor_def);
    input_data["rmi_tor"]   = Input(rmi_tor, rmi_tor_def);
    input_data["pot_tor"]   = Input(pot_tor, pot_tor_def);
    input_data["freq_tor"]  = Input(freq_tor, freq_tor_def);

    if (chem::find_section(from, key)) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else {
                Input_iter it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }
    else {
        throw Torsion_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (Cinput_iter it = input_data.begin(); it != input_data.end(); ++it) {
        if (! it->second.is_init()) {
            throw Torsion_error(it->first + " not initialized");
        }
    }

    // Check if data are sensible:

    validate();

    // Calculate reduced moment of inertia if needed:
    
    if ((rot_top.size() > 0) && rmi_tor.empty()) { 
        calc_red_imom();
    }
}

void Torsion::validate() const
{
    if (rot_axis.size() != 2) {
        throw Torsion_error("bad rot_axis size");
    }
    if (rot_axis(1) > rot.atoms.size()) {
        throw Torsion_error("bad rot_axis");
    }
    if (rot_top.size() > rot.atoms.size()) {
        throw Torsion_error("bad rot_top size");
    }
    for (arma::uword i = 0; i < rot_top.size(); ++i) {
        if (rot_top(i) > rot.atoms.size()) {
            throw Torsion_error("bad center in rot_top");
        }
    }
    if (sigma_tor.size() > 0) {
        if (arma::any(sigma_tor < 1)) {
            throw Torsion_error("bad sigma_tor");
        }
    }
    if (rmi_tor.size() > 0) {
        if (arma::any(rmi_tor <= 0.0)) {
            throw Torsion_error("bad rmi_tor");
        }
    }
    if (pot_tor.size() > 0) {
        if (arma::any(pot_tor < 0.0)) {
            throw Torsion_error("bad pot_tor");
        }
    }
    if (freq_tor.size() > 0) {
        if (arma::any(freq_tor <= 0.0)) {
            throw Torsion_error("bad freq_tor");
        }
    }
}

void Torsion::calc_red_imom()
{
    // Rotate molecule to principal axes and compute principal moments:

    rot.rotate_to_principal_axes();
    rot.principal_moments();

    // Set up axis system for rotating top:

    axis_system();
}

void Torsion::axis_system()
{
    /*
      Definition of axis system for rotating top:

            (x)       z axis is the rotating axis
            COM       x axis is perpendicular to z and goes through the
            /|          center of mass of top
           / |        y axis is perpendicular to x and z axes
        r /  | 
         /   |        C1 is the first atom center in rotating axis
        /    |        r is the vector from C1 to center of mass
       /     |_
      /      | |
      C1-------------C2 (z)
    */

    // Calculate center of mass of rotating top and work on a local copy.

    center_of_mass();
    arma::vec top_com_(top_com);

    // Set up z axis - the rotation axis - and its norm:

    z_axis = rot.xyz.row(rot_axis(1)) - rot.xyz.row(rot_axis(0));
    double z_norm = arma::norm(z_axis);

    // Find the vector from C1 to C.O.M.:

    arma::rowvec r_vec = top_com - rot.xyz.row(rot_axis(0));
    double r_norm = arma::norm(r_vec);

    // Project r vector onto z axis and find the intersection point:

    const double tol = 1.0e-12;
    double theta = arma::dot(r_vec, z_axis) / (r_norm * z_norm);
    if ((std::abs(theta) - 1.0) < tol) { // r and z are parallel
        top_com = rot.xyz.row(rot_top(0));
        r_vec = top_com - rot.xyz.row(rot_axis(0));
        r_norm = arma::norm(r_vec);
        theta = arma::dot(r_vec, z_axis) / (r_norm * z_norm);
    }
    top_origo = rot.xyz.row(rot_axis(0)) + theta * z_axis * r_norm / z_norm;

    // Set up x axis:

    x_axis = top_com - top_origo;
    double x_norm = arma::norm(x_axis);

    // Check if x and z axes are perpendicular:

    double xz_angle = arma::dot(x_axis, z_axis) / (x_norm * z_norm);
    chem::Assert(std::abs(xz_angle) < tol,
                 Torsion_error("x and z axes are not parallel"));

    // Generate y axis perpendicular to x and z axes:

    y_axis = arma::cross(z_axis, x_axis);
    double y_norm = arma::norm(y_axis);

    // Finalize by generating unit vectors:

    x_axis /= x_norm;
    y_axis /= y_norm;
    z_axis /= z_norm;
}

void Torsion::center_of_mass()
{
    double totmass = 0.0;
    for (std::size_t i = 0; i < rot.atoms.size(); ++i) {
        totmass += rot.atoms[i].atomic_mass;
    }
    for (arma::uword j = 0; j < rot.xyz.n_cols; ++j) {
        double sum = 0.0;
        for (arma::uword i = 0; i < rot_top.size(); ++i) {
            sum += rot.atoms[i].atomic_mass * rot.xyz(i,j);
        }
        top_com(j) = sum / totmass;
    }
}
