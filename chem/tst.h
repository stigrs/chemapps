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

#ifndef CHEM_TST_H
#define CHEM_TST_H

#include <chem/molecule.h>
#include <chem/thermodata.h>
#include <chem/tunnel.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// Error reporting:

struct Tst_error : std::runtime_error {
    Tst_error(std::string s) : std::runtime_error(s) {}
};

//
// Class providing Transition State Theory (TST).
//
// Note: Currently, only conventional TST is implemented. It is recommended
// to use Polyrate if variational TST is needed.
//
class Tst {
public:
    Tst(std::istream& from,
        std::ostream& to       = std::cout,
        const std::string& key = "TST",
        bool verbose           = false);

    ~Tst() {}

    // Calculate rate coefficients.
    void rate() const;

    // Calculate rate coefficients using conventional TST.
    void conventional(std::ostream& to = std::cout) const;

    // Calculate rate coefficient for the given temperature.
    double rate_coeff(double temp = 298.15) const;

private:
    // Calculate rate coefficient for the given temperature using conventional
    // TST.
    double rate_conventional(double temp = 298.15) const;

    enum Method_t { Conventional };
    enum Reaction_t { Unimolecular, Bimolecular };

    Method_t method     = Conventional;  // TST method
    Reaction_t reaction = Bimolecular;   // reaction type

    std::unique_ptr<Molecule> ra;  // reactant A
    std::unique_ptr<Molecule> rb;  // reactant B
    std::unique_ptr<Molecule> ts;  // transition state

    std::unique_ptr<Tunnel> kappa;   // tunneling correction
    std::unique_ptr<Thermodata> td;  // thermochemistry parameters

    double en_barrier;  // reaction barrier (kJ/mol)
    int rxn_sigma;      // reaction symmetry number
};

inline void Tst::rate() const
{
    switch (method) {
    case Conventional:
    default:
        conventional();
    }
}

inline double Tst::rate_coeff(double temp) const
{
    switch (method) {
    case Conventional:
    default:
        return rate_conventional(temp);
    }
}

#endif  // CHEM_TST_H
