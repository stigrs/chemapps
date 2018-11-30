// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_TST_H
#define CHEM_TST_H

#include <chem/molecule.h>
#include <chem/thermodata.h>
#include <chem/tunnel.h>
#include <iostream>

namespace Chem {

// Class providing Transition State Theory (TST).
//
// Note: Currently, only conventional TST is implemented. It is recommended
// to use Polyrate if variational TST is needed.
//
class Tst {
public:
    Tst(std::istream& from,
        std::ostream& to = std::cout,
        const std::string& key = "TST",
        bool verbose = false);

    // Calculate rate coefficients.
    void rate(std::ostream& to = std::cout) const;

    // Calculate rate coefficients using conventional TST.
    void conventional(std::ostream& to = std::cout) const;

    // Calculate rate coefficient for the given temperature.
    double rate_coeff(double temp = 298.15) const;

    // Calculate tunneling correction.
    double tunneling(double temp = 298.15) const;

private:
    // Calculate rate coefficient for the given temperature using conventional
    // TST.
    double rate_conventional(double temp = 298.15) const;

    enum Method_t { Conventional };
    enum Reaction_t { Unimolecular, Bimolecular };

    Method_t method = Conventional;    // TST method
    Reaction_t reaction = Bimolecular; // reaction type

    Thermodata td; // thermochemistry parameters
    Tunnel kappa;  // tunneling correction

    Molecule ra; // reactant A
    Molecule rb; // reactant B
    Molecule ts; // transition state

    double en_barrier; // reaction barrier (kJ/mol)
    int sigma_rxn;     // reaction symmetry number
};

inline void Tst::rate(std::ostream& to) const
{
    switch (method) {
    case Conventional:
    default:
        conventional(to);
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

inline double Tst::tunneling(double temp) const { return kappa.factor(temp); }

} // namespace Chem

#endif // CHEM_TST_H

