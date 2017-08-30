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
#include <chem/tunnel.h>
#include <chem/utils.h>
#include <map>

Tunnel::Tunnel(std::istream& from, const std::string& key)
{
    // Read input data:
    std::string method_def = "Wigner";
    std::string method_str;

    std::map<std::string, Input> input_data;
    input_data["method"]     = Input(method_str, method_def);
    input_data["freq_im"]    = Input(freq_im);
    input_data["en_barrier"] = Input(en_barrier, 0.0);
    input_data["en_rxn"]     = Input(en_rxn, 0.0);

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
        throw Tunnel_error("cannot find " + key + " section");
    }

    // Check if initialized:
    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Tunnel_error(it->first + " not initialized");
        }
    }

    // Set tunneling correction method:
    if (method_str == "Wigner") {
        method = Wigner;
    }
    else if (method_str == "Eckart") {
        method = Eckart;
    }
    else {
        throw Tunnel_error("unknown tunneling correction: " + method_str);
    }
}
