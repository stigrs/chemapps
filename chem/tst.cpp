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

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4100)  // unreferenced formal parameter
#endif                           // _MSC_VER

#include <chem/datum.h>
#include <chem/input.h>
#include <chem/thermochem.h>
#include <chem/tst.h>
#include <chem/utils.h>
#include <armadillo>
#include <gsl/gsl>
#include <map>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER

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

    // Initialize tunneling correction:
    kappa = std::make_unique<Tunnel>(from);

    // Initialize thermochemistry parameters:
    td = std::make_unique<Thermodata>(from);
}

void Tst::conventional(std::ostream& to) const
{
    arma::vec temp     = td->get_temperature();
    arma::vec pressure = arma::zeros(1);

    chem::thermochemistry(*ra, temp, pressure, false, to);
    if (reaction == Bimolecular) {
        chem::thermochemistry(*rb, temp, pressure, false, to);
    }
    chem::thermochemistry(*ts, temp, pressure, false, to);

    chem::Format<char> line;
    line.width(37).fill('=');

    to << "Conventional Transition State Theory:\n" << line('=') << "\n\n";
    if (reaction == Bimolecular) {
        to << "Reaction Rate Coefficients [cm^3 molecule^-1 s^-1]:\n";
    }
    else if (reaction == Unimolecular) {
        to << "Reaction Rate Coefficients [s^-1]:\n";
    }

    line.width(59).fill('-');
    if (kappa->get_method() == "Eckart") {
        to << line('-') << '\n'
           << "T/K\t Wigner\t Eckart  TST\t     TST/Wigner  TST/Eckart\n"
           << line('-') << '\n';
    }
    else if (kappa->get_method() == "Wigner") {
        to << line('-') << '\n'
           << "T/K\t Wigner\t TST\t     TST/Wigner\n"
           << line('-') << '\n';
    }
    else {
        to << line('-') << '\n' << "T/K\t TST\n" << line('-') << '\n';
    }

    chem::Format<double> fix7;
    fix7.fixed().width(7).precision(2);

    chem::Format<double> fix6;
    fix6.fixed().width(6).precision(2);

    chem::Format<double> sci;
    sci.scientific().width(10).precision(4);

    for (arma::uword i = 0; i < temp.size(); ++i) {
        double ktst = rate_conventional(temp(i));
        double wig  = kappa->wigner(temp(i));
        if (kappa->get_method() == "Eckart") {
            double eck = kappa->eckart(temp(i));
            to << fix7(temp(i)) << "  " << fix6(wig) << "  " << fix6(eck)
               << "  " << sci(ktst) << "  " << sci(ktst * wig) << "  "
               << sci(ktst * eck) << '\n';
        }
        else if (kappa->get_method() == "Wigner") {
            to << fix7(temp(i)) << "  " << fix6(wig) << "  " << sci(ktst)
               << "  " << sci(ktst * wig) << '\n';
        }
        else {
            to << fix7(temp(i)) << "  " << sci(ktst) << '\n';
        }
    }
    to << line('-') << '\n';
}

double Tst::rate_conventional(double temp) const
{
    Expects(temp > 0.0);

    double qts = chem::qtot(*ts, temp, 0.0, false, "V=0");
    double qa  = chem::qtot(*ra, temp, 0.0, false, "V=0");
    double qb  = 1.0;
    if (reaction == Bimolecular) {
        qb = chem::qtot(*rb, temp, 0.0, false, "V=0");
    }
    double ktst = rxn_sigma * datum::mega * datum::k * temp / datum::h;
    ktst *= (qts / (qa * qb))
            * std::exp(-en_barrier * datum::kilo / (datum::R * temp));
    return ktst;
}
