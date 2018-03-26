////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

#ifndef CHEM_MOPAC_H
#define CHEM_MOPAC_H

#include <chem/molecule.h>
#include <srs/array.h>
#include <iostream>
#include <stdexcept>
#include <string>

//-----------------------------------------------------------------------------

// Error reporting:

struct Mopac_error : std::runtime_error {
    Mopac_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

// Wrapper class for running Mopac calculations.
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
    void get_xyz(srs::dmatrix& xyz) const;

private:
    // Create Mopac input file.
    void write_dat(const Molecule& mol) const;

    // Write Cartesian coordinates in Mopac format.
    void write_xyz(std::ostream& to, const Molecule& mol) const;

    std::string version;   // Mopac version
    std::string keywords;  // list of Mopac keywords
    std::string jobname;   // Mopac job name
    int opt_geom;          // flag to specify geometry optimization
};

inline Mopac::Mopac()
{
    version  = "mopac5022mn";
    keywords = "PM6-D EF GEO-OK PRECISE";
    jobname  = "mopac";
    opt_geom = 1;  // perform geometry optimization
}

inline Mopac::Mopac(std::istream& from, const std::string& key)
{
    init(from, key);
}

#endif  // CHEM_MOPAC_H
