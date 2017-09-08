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

#ifndef CHEM_COLLISION_H
#define CHEM_COLLISION_H

#include <chem/mol_formula.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Error reporting:

struct Collision_error : std::runtime_error {
    Collision_error(std::string s) : std::runtime_error(s) {}
};

//
// Class providing methods for computing collision integrals and Lennard-Jones
// collision rates.
//
class Collision {
public:
    Collision(std::istream& from, const std::string& key = "Collision");

private:
    // Populate database with local Lennard-Jones collision diameter values.
    void set_sigma_local_values();

    // Populate database with local Lennard-Jones well depth values.
    void set_epsilon_local_values();

    enum Coll_model_t { generic, brw84, brw90a, brw90b };  // collision models
    enum Coll_omega22_t { troe, forst };  // collision integral equations

    Coll_model_t coll_model;
    Coll_omega22_t coll_integral;

    double mass_bath;     // mass of bath gas in amu
    double mass_mol;      // mass of molecule in amu
    double epsilon_bath;  // LJ well depth of bath gas in kelvin
    double epsilon_mol;   // LJ well depth of molecule in kelvin
    double sigma_bath;    // LJ collision diam. of bath gas in angstrom
    double sigma_mol;     // LJ collision diam. of molecule in angstrom
    double temperature;   // temperature in kelvin
    double number_vibr;   // number of vibrational modes of molecule
    double vibr_avg;      // average vibrational frequency of molecule in cm-1
    double vibr_high;     // highest vibrational frequency of molecule in cm-1
    double coll_energy;   // collision energy in cm-1
    double h_factor;

    std::vector<double> sigma_loc_val;     // local sigma values
    std::vector<double> epsilon_loc_val;   // local epsilon values
    std::vector<Mol_formula> mol_formula;  // molecular formula of collider
};

#endif  // CHEM_COLLISION_H
