// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/units.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

//------------------------------------------------------------------------------

enum Chem::Units::Type have;
enum Chem::Units::Type want;

void usage(const std::string& prog);
double convert(const double value);

//------------------------------------------------------------------------------

//
// Program for conversion of energies.
//
// Note: Validated against http://physics.nist.gov/cuu/Constants/energy.html.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 4) {
        usage(args[0]);
        return 1;
    }

    double value;
    std::istringstream iss(args[1]);
    iss >> value;
    if (!iss) {
        std::cerr << "Bad value: " << args[1] << '\n';
        usage(args[0]);
        return 1;
    }

    try {
        have = Chem::Units::lexer(args[2]);
        want = Chem::Units::lexer(args[3]);

        std::cout << "You have: " << value << " " << args[2] << '\n'
                  << "You want: " << convert(value) << " " << args[3] << '\n';
    }
    catch (std::exception&) {
        std::cerr << "Cannot convert " << value << " " << args[2] << " to "
                  << args[3] << '\n';
        usage(args[0]);
        return 1;
    }
}

//------------------------------------------------------------------------------

void usage(const std::string& prog)
{
    std::cerr << "Usage: " << prog << " value have_unit want_unit\n\n";
    Chem::Units::print(std::cerr);
}

double convert(const double value)
{
    namespace Pc = Numlib::Constants;

    double ans = value;
    switch (have) {
    case Chem::Units::kJ_mol:
        switch (want) {
        case Chem::Units::kcal_mol:
            ans /= Pc::cal_to_J;
            break;
        case Chem::Units::icm:
            ans /= Pc::icm_to_kJ;
            break;
        case Chem::Units::hartree:
            ans *= Pc::E_h * 1.0e-3 / Pc::N_A;
            break;
        case Chem::Units::kelvin:
            ans /= Pc::icm_to_kJ;
            ans *= Pc::icm_to_K;
            break;
        case Chem::Units::eV:
            ans /= 1.0e-3 * Pc::eV * Pc::N_A;
            break;
        case Chem::Units::kJ_mol:
            break;
        default:
            throw std::runtime_error("bad unit conversion");
        }
        break;
    case Chem::Units::kcal_mol:
        switch (want) {
        case Chem::Units::kJ_mol:
            ans *= Pc::cal_to_J;
            break;
        case Chem::Units::icm:
            ans *= Pc::cal_to_J / Pc::icm_to_kJ;
            break;
        case Chem::Units::hartree:
            ans *= Pc::cal_to_J;
            ans *= Pc::E_h * 1.0e-3 / Pc::N_A;
            break;
        case Chem::Units::kelvin:
            ans *= Pc::cal_to_J / Pc::icm_to_kJ;
            ans *= Pc::icm_to_K;
            break;
        case Chem::Units::eV:
            ans *= Pc::cal_to_J;
            ans /= 1.0e-3 * Pc::eV * Pc::N_A;
            break;
        case Chem::Units::kcal_mol:
            break;
        default:
            throw std::runtime_error("bad unit conversion");
        }
        break;
    case Chem::Units::icm:
        switch (want) {
        case Chem::Units::kJ_mol:
            ans *= Pc::icm_to_kJ;
            break;
        case Chem::Units::kcal_mol:
            ans *= Pc::icm_to_kJ / Pc::cal_to_J;
            break;
        case Chem::Units::hartree:
            ans /= Pc::au_to_icm;
            break;
        case Chem::Units::kelvin:
            ans *= Pc::icm_to_K;
            break;
        case Chem::Units::au:
            ans /= Pc::au_to_icm;
            break;
        case Chem::Units::icm:
            break;
        default:
            throw std::runtime_error("bad unit conversion");
        }
        break;
    case Chem::Units::kelvin:
        switch (want) {
        case Chem::Units::icm:
            ans /= Pc::icm_to_K;
            break;
        case Chem::Units::au:
            ans /= Pc::au_to_K;
            break;
        case Chem::Units::kelvin:
            break;
        default:
            throw std::runtime_error("bad unit conversion");
        }
        break;
    case Chem::Units::au:
        switch (want) {
        case Chem::Units::icm:
            ans *= Pc::au_to_icm;
            break;
        case Chem::Units::kelvin:
            ans *= Pc::au_to_K;
            break;
        case Chem::Units::kg:
            ans *= Pc::au_to_kg;
            break;
        case Chem::Units::eV:
            ans *= Pc::E_h / Pc::eV;
            break;
        case Chem::Units::au:
            break;
        default:
            throw std::runtime_error("bad unit conversion");
        }
        break;
    case Chem::Units::eV:
        switch (want) {
        case Chem::Units::kJ_mol:
            ans *= 1.0e-3 * Pc::eV * Pc::N_A;
            break;
        case Chem::Units::kcal_mol:
            ans *= 1.0e-3 * Pc::eV * Pc::N_A;
            ans /= Pc::cal_to_J;
            break;
        case Chem::Units::hartree:
            ans *= Pc::eV / Pc::E_h;
            break;
        case Chem::Units::eV:
            break;
        default:
            throw std::runtime_error("bad unit conversion");
        }
        break;
    default:
        throw std::runtime_error("unknown unit");
    }
    return ans;
}

