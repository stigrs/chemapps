// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_MOPAC_H
#define CHEM_MOPAC_H

#include <chem/molecule.h>
#include <numlib/matrix.h>
#include <iostream>
#include <string>

namespace Chem {

// Wrapper class for running Mopac calculations.
//
class Mopac {
public:
    Mopac();

    Mopac(std::istream& from, const std::string& key = "Mopac");

    // Initialize Mopac calculation.
    void init(std::istream& from, const std::string& key = "Mopac");

    // Run Mopac calculation.
    void run(Molecule& mol) const;

    // Check SCF convergence.
    bool check_convergence() const;

    // Get heat of formation in kJ/mol from Mopac output file.
    double get_heat_of_formation() const;

    // Get optimized Cartesian coordinates.
    void get_xyz(Numlib::Mat<double>& xyz) const;

private:
    // Create Mopac input file.
    void write_dat(const Molecule& mol) const;

    // Write Cartesian coordinates in Mopac format.
    void write_xyz(std::ostream& to, const Molecule& mol) const;

    std::string version;  // Mopac version
    std::string keywords; // list of Mopac keywords
    std::string jobname;  // Mopac job name
    int opt_geom;         // flag to specify geometry optimization
};

inline Mopac::Mopac()
{
    version = "mopac5022mn";
    keywords = "PM6-D EF GEO-OK PRECISE";
    jobname = "mopac";
    opt_geom = 1; // perform geometry optimization
}

inline Mopac::Mopac(std::istream& from, const std::string& key)
{
    init(from, key);
}

} // namespace Chem

#endif // CHEM_MOPAC_H

