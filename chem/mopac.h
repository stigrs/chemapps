/**
   @file mopac.h
   
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

#ifndef CHEM_MOPAC_H
#define CHEM_MOPAC_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <chem/molecule.h>

//-----------------------------------------------------------------------------

// Error reporting:

struct Mopac_error : std::runtime_error {
    Mopac_error(std::string s) : std::runtime_error(s) { }
};

//-----------------------------------------------------------------------------

/// Wrapper class for running Mopac calculations.
class Mopac {
public:
    Mopac();

    Mopac(std::istream& from, const std::string& key);

    ~Mopac() { }

    /// Run Mopac calculation.
    void run(Molecule& mol) const;
    
private:
    /// Create Mopac input file.
    void write_dat(const Molecule& mol, const std::string& dat_file) const;

    /// Write Cartesian coordinates in Mopac format.
    void write_xyz(std::ostream& to, const Molecule& mol) const;
    
    std::string version;  ///< Mopac version
    std::string keywords; ///< list of Mopac keywords
    std::string jobname;  ///< Mopac job name
    int opt_geom;         ///< flag to specify geometry optimization
};

inline Mopac::Mopac()
{
    version = "mopac5022mn";
    keywords = "PM6-D GEO-OK PRECISE";
    jobname = "mopac";
    opt_geom = 1; // true
}

#endif /* CHEM_MOPAC_H */

