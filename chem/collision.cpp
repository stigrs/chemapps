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
    // input_data["mol_formula"]   = Input(mol_formula, mol_formula);
    input_data["h_factor"] = Input(h_factor, 1.0);

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