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

#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

struct Cambi_error : std::runtime_error {
    Cambi_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Structure for holding input data.
struct VdW_data {
    std::string name; // name of species
    std::string type; // type of species (atom or molecule)
    double alpha;     // polarizability (angstrom**3)
    int n_el_int;     // number of inner electrons
    int n_el_ext;     // number of outer electrons
};

//------------------------------------------------------------------------------

// Forward declarations:

double n_el_eff(int species);
double rm_eq_dist();
double c6_eff();
double lj_well_depth();

void init();
void read_data(std::istream& is, const std::streamoff& pos, VdW_data& data);
void print_data(const VdW_data& data);

//------------------------------------------------------------------------------

// Global declarations:

std::ifstream from;
VdW_data a_data;
VdW_data b_data;

//------------------------------------------------------------------------------

// Program for calculation of van der Waals interaction potential parameters
// using the generalized correlations presented in the paper:
//
// Cambi, R.; Cappelletti, D.; Liuti, G.; Pirani, F. J. Chem. Phys.,
// 1991, vol. 95, pp. 1852-1861.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Calculate van der Waals interaction");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>());
    // clang-format on

    auto result = options.parse(argc, argv);

    std::string input_file;

    if (result.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (result.count("file")) {
        input_file = result["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        Stdutils::fopen(from, input_file);

        init();

        const double meV_to_K = 11.6045;
        const double meV_to_cm = 8.065545;

        double sigma = rm_eq_dist() / std::pow(2.0, 1. / 6.);
        double epsilon = lj_well_depth();

        std::cout << "van der Waals interaction potential parameters: \n"
                  << "  Neff_A:  " << n_el_eff(1) << '\n'
                  << "  Neff_B:  " << n_el_eff(2) << '\n'
                  << "  R_m:     " << rm_eq_dist() << " angstrom\n"
                  << "  sigma:   " << sigma << " angstrom\n"
                  << "  C6_eff:  " << c6_eff() << " meV angstrom**6\n"
                  << "  epsilon: " << epsilon << " meV\n"
                  << "           " << epsilon * meV_to_K << " K\n"
                  << "           " << epsilon * meV_to_cm << " cm**-1\n";
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

double n_el_eff(int species)
{
    double n_tot;
    double n_eff;
    double n_int;
    double n_ext;
    std::string type;
    switch (species) {
    case 1:
        n_int = a_data.n_el_int;
        n_ext = a_data.n_el_ext;
        type = a_data.type;
        break;
    case 2:
        n_int = b_data.n_el_int;
        n_ext = b_data.n_el_ext;
        type = b_data.type;
        break;
    default:
        throw Cambi_error("unknown species number: " + std::to_string(species));
    }
    n_tot = n_int + n_ext;

    if (type == "atom") {
        n_eff = 1.0 + (1.0 - n_ext / n_int) * std::pow(n_int / n_tot, 2.0);
        n_eff = n_eff * n_ext;
    }
    else { // molecule
        n_eff = 1.0 - n_int * n_ext / std::pow(n_tot, 2.0);
        n_eff = n_eff * n_tot;
    }
    return n_eff;
}

double rm_eq_dist()
{
    const double gamma = 0.095;
    const double rm_coeff = 1.767;

    return rm_coeff *
           (std::pow(a_data.alpha, 1. / 3.) + std::pow(b_data.alpha, 1. / 3.)) /
           std::pow(a_data.alpha * b_data.alpha, gamma);
}

double c6_eff()
{
    double alphaA_N = std::sqrt(a_data.alpha / n_el_eff(1));
    double alphaB_N = std::sqrt(b_data.alpha / n_el_eff(2));

    return 15.7e+3 * a_data.alpha * b_data.alpha / (alphaA_N + alphaB_N);
}

double lj_well_depth()
{
    return 0.720 * c6_eff() / std::pow(rm_eq_dist(), 6.0);
}

void init()
{
    auto pos = Stdutils::find_token(from, "SpeciesA");
    if (pos != -1) {
        try {
            read_data(from, pos, a_data);
        }
        catch (std::exception& e) {
            std::cerr << e.what() << '\n';
            throw Cambi_error("failed to read data on species A");
        }
    }
    else {
        throw Cambi_error("could not find section 'SpeciesA'");
    }

    pos = Stdutils::find_token(from, "SpeciesB");
    if (pos != -1) {
        try {
            read_data(from, pos, b_data);
        }
        catch (std::exception& e) {
            std::cerr << e.what() << '\n';
            throw Cambi_error("failed to read data on species B");
        }
    }
    else {
        throw Cambi_error("could not find section 'SpeciesB'");
    }
    print_data(a_data);
    print_data(b_data);
}

void read_data(std::istream& is, const std::streamoff& pos, VdW_data& data)
{
    using namespace Stdutils;

    get_token_value(is, pos, "name", data.name);
    get_token_value(is, pos, "type", data.type);
    get_token_value(is, pos, "polarizability", data.alpha);
    get_token_value(is, pos, "n_el_int", data.n_el_int);
    get_token_value(is, pos, "n_el_ext", data.n_el_ext);

    // Check if data are sensible:

    if ((data.type != "atom") && (data.type != "molecule")) {
        throw Cambi_error("type has bad value: " + data.type);
    }
    if (data.alpha < 0) {
        throw Cambi_error("polarizability has bad value: " +
                          std::to_string(data.alpha));
    }
    if (data.n_el_int < 0) {
        throw Cambi_error("n_el_int has bad value: " +
                          std::to_string(data.n_el_int));
    }
    if (data.n_el_ext < 0) {
        throw Cambi_error("n_el_ext has bad value: " +
                          std::to_string(data.n_el_ext));
    }
}

void print_data(const VdW_data& data)
{
    std::cout << "Input data on " << data.name << ":\n"
              << "  polarizability: " << data.alpha << " angstrom**3\n";
    if (data.type == "atom") {
        std::cout << "  number of inner electrons: " << data.n_el_int << '\n'
                  << "  number of outer electrons: " << data.n_el_ext << "\n\n";
    }
    else { // molecule
        std::cout << "  number of bonding electrons: " << data.n_el_int << '\n'
                  << "  number of non-bonding electrons: " << data.n_el_ext
                  << "\n\n";
    }
}

