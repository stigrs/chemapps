/**
   @file mopac.cpp
   
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
#include <chem/mopac.h>
#include <chem/input.h>
#include <chem/utils.h>


Mopac::Mopac(std::istream& from, const std::string& key)
{
    typedef std::map<std::string, Input>::iterator       Input_iter;
    typedef std::map<std::string, Input>::const_iterator Cinput_iter;

    // Read input data:

    std::string version_def = "mopac5022mn";
    std::string keywords_def = "PM6-D GEO-OK PRECISE";
    std::string jobname_def = "mopac";
    int opt_geom_def = 1;
    
    std::map<std::string, Input> input_data;
    input_data["version"]  = Input(version, version_def);
    input_data["jobname"]  = Input(jobname, jobname_def);
    input_data["opt_geom"] = Input(opt_geom, opt_geom_def);

    if (chem::find_section(from, key)) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "keywords") {
                std::string line;
                std::getline(from, line);
                keywords = line;
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
        throw Mopac_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (Cinput_iter it = input_data.begin(); it != input_data.end(); ++it) {
        if (! it->second.is_init()) {
            throw Mopac_error(it->first + " not initialized");
        }
    }
    for (Cinput_iter it = input_data.begin(); it != input_data.end(); ++it) {
        std::cout << it->first << ": " << it->second << '\n';
    }    
}
