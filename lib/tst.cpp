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

#include <chem/thermochem.h>
#include <chem/tst.h>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <string>
#include <stdexcept>

Chem::Tst::Tst(std::istream& from,
               std::ostream& to,
               const std::string& key,
               bool verbose)
    : td(from), kappa(from)
{
    using namespace Stdutils;

    // Read input data:
    std::string method_str = "conventional";
    std::string reaction_str = "bimolecular";

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "method", method_str, method_str);
        get_token_value(from, pos, "reaction", reaction_str, reaction_str);
        get_token_value(from, pos, "en_barrier", en_barrier, 0.0);
        get_token_value(from, pos, "sigma_rxn", sigma_rxn, 1);
    }
    Assert::dynamic(en_barrier > 0.0, "bad energy barrier");
    Assert::dynamic(sigma_rxn >= 1, "bad reaction multiplicity");

    // Set TST method:
    if (method_str == "Conventional") {
        method = Conventional;
    }
    else {
        throw std::runtime_error("unknown TST method: " + method_str);
    }

    // Set reaction type:
    if (reaction_str == "Unimolecular") {
        reaction = Unimolecular;
    }
    else if (reaction_str == "Bimolecular") {
        reaction = Bimolecular;
    }
    else {
        throw std::runtime_error("unknown reaction: " + reaction_str);
    }

    // Initialize reactants and transition state:

    ra = std::make_unique<Chem::Molecule>(from, to, "ReactantA", verbose);

    if (reaction == Bimolecular) {
        rb = std::make_unique<Chem::Molecule>(from, to, "ReactantB", verbose);
    }
    ts = std::make_unique<Chem::Molecule>(from, to, "TransitionState", verbose);
}

void Chem::Tst::conventional(std::ostream& to) const
{
    Numlib::Vec<double> temp = td.get_temperature();
    Numlib::Vec<double> pressure = {0.0};

    Chem::thermochemistry(*ra, temp, pressure, false, to);
    if (reaction == Bimolecular) {
        Chem::thermochemistry(*rb, temp, pressure, false, to);
    }
    Chem::thermochemistry(*ts, temp, pressure, false, to);

    Stdutils::Format<char> line;
    line.width(37).fill('=');

    to << "Conventional Transition State Theory:\n" << line('=') << "\n\n";
    if (reaction == Bimolecular) {
        to << "Reaction Rate Coefficients [cm^3 molecule^-1 s^-1]:\n";
    }
    else if (reaction == Unimolecular) {
        to << "Reaction Rate Coefficients [s^-1]:\n";
    }

    line.width(59).fill('-');
    if (kappa.get_method() == "Eckart") {
        to << line('-') << '\n'
           << "T/K\t Wigner\t Eckart  TST\t     TST/Wigner  TST/Eckart\n"
           << line('-') << '\n';
    }
    else if (kappa.get_method() == "Wigner") {
        to << line('-') << '\n'
           << "T/K\t Wigner\t TST\t     TST/Wigner\n"
           << line('-') << '\n';
    }
    else {
        to << line('-') << '\n' << "T/K\t TST\n" << line('-') << '\n';
    }

    Stdutils::Format<double> fix7;
    fix7.fixed().width(7).precision(2);

    Stdutils::Format<double> fix6;
    fix6.fixed().width(6).precision(2);

    Stdutils::Format<double> sci;
    sci.scientific().width(10).precision(4);

    for (Index i = 0; i < temp.size(); ++i) {
        double ktst = rate_conventional(temp(i));
        double wig = kappa.wigner(temp(i));
        if (kappa.get_method() == "Eckart") {
            double eck = kappa.eckart(temp(i));
            to << fix7(temp(i)) << "  " << fix6(wig) << "  " << fix6(eck)
               << "  " << sci(ktst) << "  " << sci(ktst * wig) << "  "
               << sci(ktst * eck) << '\n';
        }
        else if (kappa.get_method() == "Wigner") {
            to << fix7(temp(i)) << "  " << fix6(wig) << "  " << sci(ktst)
               << "  " << sci(ktst * wig) << '\n';
        }
        else {
            to << fix7(temp(i)) << "  " << sci(ktst) << '\n';
        }
    }
    to << line('-') << '\n';
}

double Chem::Tst::rate_conventional(double temp) const
{
    using namespace Numlib::Constants;

    Assert::dynamic<Assert::level(2)>(temp > 0.0, "bad temperature");

    double qts = Chem::qtot(*ts, temp, 0.0, false, "V=0");
    double qa = Chem::qtot(*ra, temp, 0.0, false, "V=0");
    double qb = 1.0;
    if (reaction == Bimolecular) {
        qb = Chem::qtot(*rb, temp, 0.0, false, "V=0");
    }
    double ktst = sigma_rxn * mega * k * temp / h;
    ktst *= (qts / (qa * qb)) * std::exp(-en_barrier * kilo / (R * temp));
    return ktst;
}

