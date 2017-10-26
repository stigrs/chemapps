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

#include <chem/collision.h>
#include <chem/input.h>
#include <chem/ptable.h>
#include <chem/utils.h>
#include <gsl/gsl>
#include <map>

Collision::Collision(std::istream& from, const std::string& key)
{
    // Populate local sigma and epsilon values:
    set_sigma_local_values();
    set_epsilon_local_values();

    // Read input data:
    std::string coll_model_str;
    std::string coll_integral_str;

    std::map<std::string, Input> input_data;
    input_data["coll_model"]    = Input(coll_model_str, "generic");
    input_data["coll_integral"] = Input(coll_integral_str, "forst");
    input_data["mass_bath"]     = Input(mass_bath);
    input_data["mass_mol"]      = Input(mass_mol);
    input_data["epsilon_bath"]  = Input(epsilon_bath);
    input_data["epsilon_mol"]   = Input(epsilon_mol);
    input_data["sigma_bath"]    = Input(sigma_bath);
    input_data["sigma_mol"]     = Input(sigma_mol);
    input_data["number_vibr"]   = Input(number_vibr, 0);
    input_data["vibr_avg"]      = Input(vibr_avg, 0.0);
    input_data["vibr_high"]     = Input(vibr_high, 0.0);
    input_data["temperature"]   = Input(temperature);
    input_data["coll_energy"]   = Input(coll_energy, 0.0);
    input_data["mol_formula"]   = Input(mol_formula, mol_formula);
    input_data["h_factor"]      = Input(h_factor, 1.0);

    if (chem::find_section(from, key)) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }
    else {
        throw Collision_error("cannot find " + key + " section");
    }

    // Check if initialized:
    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Collision_error(it->first + " not initialized");
        }
    }

    // Validate input data:
    if (coll_model_str == "generic") {
        coll_model = generic;
    }
    else if (coll_model_str == "brw84") {
        coll_model = brw84;
    }
    else if (coll_model_str == "brw90a") {
        coll_model = brw90a;
    }
    else if (coll_model_str == "brw90b") {
        coll_model = brw90b;
    }
    else {
        throw Collision_error("bad coll_model: " + coll_model_str);
    }
    if (coll_integral_str == "troe") {
        coll_integral = troe;
    }
    else if (coll_integral_str == "forst") {
        coll_integral = forst;
    }
    else {
        throw Collision_error("bad coll_integral: " + coll_integral_str);
    }
    Expects(mass_bath > 0.0);
    Expects(mass_mol > 0.0);
    Expects(epsilon_bath > 0.0);
    Expects(epsilon_mol > 0.0);
    Expects(sigma_bath > 0.0);
    Expects(sigma_mol > 0.0);
    Expects(temperature > 0.0);
    if (coll_model != generic) {
        Expects(number_vibr >= 1);
    }
    if (coll_model == brw84) {
        Expects(vibr_avg > 0.1);
        Expects((h_factor > 0.0) && (h_factor <= 1.0));
    }
    if ((coll_model == brw84) || (coll_model == brw90a)) {
        Expects(coll_energy > 0.1);
    }
    if ((coll_model == brw90a) || (coll_model == brw90b)) {
        Expects(vibr_high > 0.1);
        Expects(!mol_formula.empty());
    }
}

double Collision::average_mass() const
{
    double mavg = 0.0;

    int natoms = 0;
    for (std::size_t i = 0; i < mol_formula.size(); ++i) {
        natoms += mol_formula[i].stoich;
        mavg += mol_formula[i].stoich
                * ptable::get_atomic_mass(mol_formula[i].atom);
    }
    mavg /= gsl::narrow_cast<double>(natoms);  // avr. mass of atoms in molecule

    if (coll_model == brw90a) {  // eq. 35a in Kim and Gilbert (1990)
        mavg = 1.0 / (1.0 / reduced_mass() + 1.0 / mavg);
    }
    else if (coll_model == brw90b) {  // eq. 35b in Kim and Gilbert (1990)
        mavg = 1.0 / ((1.0 / (mavg * natoms - mavg)) + (1.0 / mavg));
    }
    return mavg;
}

double Collision::sigma_local() const
{
    using namespace ptable;

    double sloc = 0.0;

    if ((coll_model == brw90a) || (coll_model == brw90b)) {
        int natoms = 0;
        for (std::size_t i = 0; i < mol_formula.size(); ++i) {
            natoms += mol_formula[i].stoich;
            sloc += mol_formula[i].stoich
                    * sigma_loc_val[get_atomic_number(mol_formula[i].atom)];
        }
        sloc /= gsl::narrow_cast<double>(natoms);
        sloc = 0.5 * (sloc + sigma_bath);
    }
    return sloc;
}

double Collision::epsilon_local() const
{
    using namespace ptable;

    double eloc = 0.0;

    if ((coll_model == brw90a) || (coll_model == brw90b)) {
        int natoms = 0;
        for (std::size_t i = 0; i < mol_formula.size(); ++i) {
            natoms += mol_formula[i].stoich;
            eloc += mol_formula[i].stoich
                    * epsilon_loc_val[get_atomic_number(mol_formula[i].atom)];
        }
        eloc /= gsl::narrow_cast<double>(natoms);
        eloc = std::sqrt(eloc * epsilon_bath);
    }
    return eloc;
}

void Collision::set_sigma_local_values()
{
    sigma_loc_val.resize(ptable::get_max_atomic_number(), 0.0);
    sigma_loc_val[ptable::get_atomic_number("H")]  = 3.0;
    sigma_loc_val[ptable::get_atomic_number("C")]  = 3.2;
    sigma_loc_val[ptable::get_atomic_number("N")]  = 3.2;
    sigma_loc_val[ptable::get_atomic_number("O")]  = 3.2;
    sigma_loc_val[ptable::get_atomic_number("S")]  = 3.4;
    sigma_loc_val[ptable::get_atomic_number("F")]  = 3.2;
    sigma_loc_val[ptable::get_atomic_number("Cl")] = 3.4;
    sigma_loc_val[ptable::get_atomic_number("Br")] = 3.6;
    sigma_loc_val[ptable::get_atomic_number("I")]  = 4.0;
}

void Collision::set_epsilon_local_values()
{
    epsilon_loc_val.resize(ptable::get_max_atomic_number(), 0.0);
    epsilon_loc_val[ptable::get_atomic_number("H")]  = 6.5;
    epsilon_loc_val[ptable::get_atomic_number("C")]  = 20.3;
    epsilon_loc_val[ptable::get_atomic_number("N")]  = 20.3;
    epsilon_loc_val[ptable::get_atomic_number("O")]  = 20.3;
    epsilon_loc_val[ptable::get_atomic_number("S")]  = 120.0;
    epsilon_loc_val[ptable::get_atomic_number("F")]  = 20.3;
    epsilon_loc_val[ptable::get_atomic_number("Cl")] = 120.0;
    epsilon_loc_val[ptable::get_atomic_number("Br")] = 190.0;
    epsilon_loc_val[ptable::get_atomic_number("I")]  = 230.0;
}