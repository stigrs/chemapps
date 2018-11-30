// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_GAUSSIAN_H
#define CHEM_GAUSSIAN_H

#include <chem/molecule.h>
#include <iostream>
#include <string>

namespace Chem {

// Wrapper class for running Gaussian calculations.
//
class Gaussian {
public:
    Gaussian();

    Gaussian(std::istream& from, const std::string& key = "Gaussian");

    // Initialize Gaussian calculation.
    void init(std::istream& from, const std::string& key = "Gaussian");

    // Run Gaussian calculation.
    void run(Molecule& mol) const;

private:
    // Create Gaussian input file.
    void write_com(const Molecule& mol) const;

    std::string version;  // Gaussian version
    std::string keywords; // list of Gaussian keywords
    std::string jobname;  // Gaussian job name
    int nprocshared;      // number of processors
    bool nosave;          // flag to specify if chk file should be saved
};

inline Gaussian::Gaussian()
{
    version = "rung09";
    keywords = "opt freq hf/sto-3g";
    jobname = "gauss";
    nprocshared = 1;
    nosave = true;
}

inline Gaussian::Gaussian(std::istream& from, const std::string& key)
{
    init(from, key);
}

} // namespace Chem

#endif // CHEM_GAUSSIAN_H

