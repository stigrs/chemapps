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

#include <chem/input.h>
#include <chem/tst.h>
#include <chem/utils.h>
#include <map>

Tst::Tst(std::istream& from,
         std::ostream& to,
         const std::string& key,
         bool verbose)
{
    // Read input data:
    std::string method_def   = "conventional";
    std::string reaction_def = "bimolecular";
    std::string method_str;
    std::string reaction_str;

    std::map<std::string, Input> input_data;
    input_data["method"]     = Input(method_str, method_def);
    input_data["reaction"]   = Input(reaction_str, reaction_def);
    input_data["en_barrier"] = Input(en_barrier);
    input_data["rxn_sigma"]  = Input(rxn_sigma, 1);

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
        throw Tst_error("cannot find " + key + " section");
    }

    // Check if initialized:
    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Tst_error(it->first + " not initialized");
        }
    }

    // Set TST method:
    if (method_str == "Conventional") {
        method = Conventional;
    }
    else {
        throw Tst_error("unknown TST method: " + method_str);
    }

    // Set reaction type:
    if (reaction_str == "Unimolecular") {
        reaction = Unimolecular;
    }
    else if (reaction_str == "Bimolecular") {
        reaction = Bimolecular;
    }
    else {
        throw Tst_error("unknown reaction: " + reaction_str);
    }

    // Initialize reactants and transition state:
    ra = std::make_unique<Molecule>(from, to, "ReactantA", verbose);
    if (reaction == Bimolecular) {
        rb = std::make_unique<Molecule>(from, to, "ReactantB", verbose);
    }
    ts = std::make_unique<Molecule>(from, to, "TransitionState", verbose);
}
