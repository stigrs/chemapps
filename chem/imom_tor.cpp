/**
   @file imom_tor.cpp
   
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
#include <chem/imom_tor.h>
#include <chem/input.h>
#include <chem/utils.h>


void Imom_tor::init(std::istream& from, const std::string& key)
{
    typedef std::map<std::string, Input>::iterator       Input_iter;
    typedef std::map<std::string, Input>::const_iterator Cinput_iter;

    // Read input data:

    arma::uvec rot_axis_def = arma::zeros<arma::uvec>(2);
    arma::uvec rot_top_def = arma::zeros<arma::uvec>(0);

    std::map<std::string, Input> input_data;
    input_data["rot_axis"] = Input(rot_axis, rot_axis_def);
    input_data["rot_top"]  = Input(rot_top, rot_top_def);

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
        throw Imom_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (Cinput_iter it = input_data.begin(); it != input_data.end(); ++it) {
        if (! it->second.is_init()) {
            throw Imom_error(it->first + " not initialized");
        }
    }

    // Check if data are sensible:

    if (rot_axis.size() != 2) {
        throw Imom_error("bad rot_axis size");
    }
    if (rot_axis(1) > rot.atoms.size()) {
        throw Imom_error("bad rot_axis");
    }
    if (rot_top.size() > rot.atoms.size()) {
        throw Imom_error("bad rot_top size");
    }
    for (arma::uword i = 0; i < rot_top.size(); ++i) {
        if (rot_top(i) > rot.atoms.size()) {
            throw Imom_error("bad center in rot_top");
        }
    }
}
