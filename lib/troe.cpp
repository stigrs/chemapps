////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013-2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/troe.h>
#include <srs/utils.h>
#include <map>

Troe::Troe(std::istream& from, Molecule& mol_) : mol(mol_)
{
    int pot_type_tmp;

    std::map<std::string, srs::Input> input_data;
    input_data["pot_type"]    = srs::Input(pot_type_tmp, 1);
    input_data["e_barrier"]   = srs::Input(e_barrier);
    input_data["imom_ratio"]  = srs::Input(imom_ratio, 1.0);
    input_data["n_morse_osc"] = srs::Input(n_morse_osc, 0);

    if (srs::find_section(from, "Troe")) {
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
        throw Troe_error("could not find Troe section");
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Troe_error(it->first + " not initialized");
        }
    }

    // Check if data are sensible:

    if (pot_type_tmp == 1) {
        pot_type = type1;
    }
    else if (pot_type_tmp == 2) {
        pot_type = type2;
    }
    else {
        throw Troe_error("bad potential type");
    }
    if (e_barrier <= 0.0) {
        throw Troe_error("bad energy barrier");
    }
    if (imom_ratio < 0.0) {
        throw Troe_error("bad moment of inertia ratio");
    }
    if (n_morse_osc < 0) {
        throw Troe_error("bad number of Morse oscillators");
    }
}