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

double Tunnel::eckart(double temp) const
{
    // The implementation is based on the following papers:
    //  1. Eckart, E. Phys. Rev., 1962, vol. 35, p. 1303.
    //  2. Brown, R. L. J. Research NIST, 1981, vol. 86, p. 357.
    //  3. Johnston, H. S.; Heicklen, J. J. Phys. Chem., 1962, vol. 66, p. 532.

    // Algorithm:
    // ----------
    // Brown (1981) introduced a new variable
    //
    //     epsilon = (energy - en_barrier) / (boltzmann * temp)
    //
    // in order to evaluate the integral yielding the tunneling correction.
    // When epsilon gets large, the transmission probability, kappa, approaches
    // unity. This method uses a Gaussian formula for the part of the integral
    // where kappa < 1. The remainder where kappa is ca. 1 is evaluated
    // analytically. The energy at which this happens is called epsilon_b and
    // kappa is called kappa_b. However, epsilon_b is kept below a certain
    // constant value epsilon_max. The parameters kappa_b and epsilon_max have
    // been adjusted to give best possible agreement with the values presented
    // in Table I in the paper of Johnston and Heicklen (1962). The difference
    // is less than 1 percent and therefore should give better results than
    // the program of Brown (1962). This has not been fully tested though.

    Expects(temp > 0.0);
    return temp;
}