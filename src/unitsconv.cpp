////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009-2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/units.h>
#include <srs/datum.h>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>


//------------------------------------------------------------------------------

enum Units::Type have;
enum Units::Type want;

void usage(const char* prog);
double convert(const double value);

//------------------------------------------------------------------------------

//
// Program for conversion of energies.
//
// Note: Validated against http://physics.nist.gov/cuu/Constants/energy.html.
//
int main(int argc, char* argv[])
{
    if (argc != 4) {
        usage(argv[0]);
        return 1;
    }

    double value;
    std::istringstream iss(argv[1]);
    iss >> value;
    if (!iss) {
        std::cerr << "Bad value: " << argv[1] << '\n';
        usage(argv[0]);
        return 1;
    }

    try {
        have = Units::lexer(argv[2]);
        want = Units::lexer(argv[3]);

        std::cout << "You have: " << value << " " << argv[2] << '\n'
                  << "You want: " << convert(value) << " " << argv[3] << '\n';
    }
    catch (std::exception&) {
        std::cerr << "Cannot convert " << value << " " << argv[2] << " to "
                  << argv[3] << '\n';
        usage(argv[0]);
        return 1;
    }
}

//------------------------------------------------------------------------------

void usage(const char* prog)
{
    std::cerr << "Usage: " << prog << " value have_unit want_unit\n\n";
    Units::print(std::cerr);
}

double convert(const double value)
{
    double ans = value;
    switch (have) {
    case Units::kJ_mol:
        switch (want) {
        case Units::kcal_mol:
            ans /= datum::cal_to_J;
            break;
        case Units::icm:
            ans /= datum::icm_to_kJ;
            break;
        case Units::hartree:
            ans *= datum::E_h * 1.0e-3 / datum::N_A;
            break;
        case Units::kelvin:
            ans /= datum::icm_to_kJ;
            ans *= datum::icm_to_K;
            break;
        case Units::eV:
            ans /= 1.0e-3 * datum::eV * datum::N_A;
            break;
        case Units::kJ_mol:
            break;
        default:
            throw Unit_error("bad unit conversion");
        }
        break;
    case Units::kcal_mol:
        switch (want) {
        case Units::kJ_mol:
            ans *= datum::cal_to_J;
            break;
        case Units::icm:
            ans *= datum::cal_to_J / datum::icm_to_kJ;
            break;
        case Units::hartree:
            ans *= datum::cal_to_J;
            ans *= datum::E_h * 1.0e-3 / datum::N_A;
            break;
        case Units::kelvin:
            ans *= datum::cal_to_J / datum::icm_to_kJ;
            ans *= datum::icm_to_K;
            break;
        case Units::eV:
            ans *= datum::cal_to_J;
            ans /= 1.0e-3 * datum::eV * datum::N_A;
            break;
        case Units::kcal_mol:
            break;
        default:
            throw Unit_error("bad unit conversion");
        }
        break;
    case Units::icm:
        switch (want) {
        case Units::kJ_mol:
            ans *= datum::icm_to_kJ;
            break;
        case Units::kcal_mol:
            ans *= datum::icm_to_kJ / datum::cal_to_J;
            break;
        case Units::hartree:
            ans /= datum::au_to_icm;
            break;
        case Units::kelvin:
            ans *= datum::icm_to_K;
            break;
        case Units::au:
            ans /= datum::au_to_icm;
            break;
        case Units::icm:
            break;
        default:
            throw Unit_error("bad unit conversion");
        }
        break;
    case Units::kelvin:
        switch (want) {
        case Units::icm:
            ans /= datum::icm_to_K;
            break;
        case Units::au:
            ans /= datum::au_to_K;
            break;
        case Units::kelvin:
            break;
        default:
            throw Unit_error("bad unit conversion");
        }
        break;
    case Units::au:
        switch (want) {
        case Units::icm:
            ans *= datum::au_to_icm;
            break;
        case Units::kelvin:
            ans *= datum::au_to_K;
            break;
        case Units::kg:
            ans *= datum::au_to_kg;
            break;
        case Units::eV:
            ans *= datum::E_h / datum::eV;
            break;
        case Units::au:
            break;
        default:
            throw Unit_error("bad unit conversion");
        }
        break;
    case Units::eV:
        switch (want) {
        case Units::kJ_mol:
            ans *= 1.0e-3 * datum::eV * datum::N_A;
            break;
        case Units::kcal_mol:
            ans *= 1.0e-3 * datum::eV * datum::N_A;
            ans /= datum::cal_to_J;
            break;
        case Units::hartree:
            ans *= datum::eV / datum::E_h;
            break;
        case Units::eV:
            break;
        default:
            throw Unit_error("bad unit conversion");
        }
        break;
    default:
        throw Unit_error("unknown unit");
    }
    return ans;
}
