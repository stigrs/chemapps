///////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2004-2018 Stig Rune Sellevag. All rights reserved.
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

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4505)  // caused by boost/program_options.hpp>
#endif                           // _MSC_VER

#include <chem/gauss_data.h>
#include <srs/array.h>
#include <srs/datum.h>
#include <srs/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER


//------------------------------------------------------------------------------

// Global variables:

int natoms;

double energy0;
double s_corr;

bool change_sign;
bool reverse_mep;
bool get_hess;

srs::dvector mep;
srs::dvector geom;
srs::dvector grad;
srs::dvector hess;

//------------------------------------------------------------------------------

struct Error : std::runtime_error {
    Error(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

void create_fu31(const Gauss_filetype filetype);
Gauss_filetype get_filetype(const char* filename);

//------------------------------------------------------------------------------

//
// Prepares output from a Gaussian 94/98/03 IRC calculation for
// input to Polyrate's unit31.
//
int main(int argc, char* argv[])
{
    // Set default values:

    energy0 = 0.0;
    s_corr  = 0.0;

    change_sign = false;
    reverse_mep = false;
    get_hess    = false;

    // Parse command line arguments:

    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file")
        ("energy,e", po::value<double>(), "energy of reactants (0.0)")
        ("corr,c", po::value<double>(), "correction to SMEP (0.0)")
        ("sign,s", po::bool_switch()->default_value(false), "change sign of SMEP values (false)")
        ("reverse,r", po::bool_switch()->default_value(false), "reverse MEP values (false")
        ("hess,h", po::bool_switch()->default_value(false), "get Hessians (false");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    std::string input_file;

    if (vm.find("help") != vm.end()) {
        std::cout << options << '\n';
        return 0;
    }
    if (vm.find("file") != vm.end()) {
        input_file = vm["file"].as<std::string>();
    }
    else {
        std::cerr << options << '\n';
        return 1;
    }
    if (vm.find("energy") != vm.end()) {
        energy0 = vm["energy"].as<double>();
    }
    if (vm.find("corr") != vm.end()) {
        s_corr = vm["corr"].as<double>();
    }
    if (vm.find("sign") != vm.end()) {
        change_sign = vm["sign"].as<bool>();
    }
    if (vm.find("reverse") != vm.end()) {
        reverse_mep = vm["reverse"].as<bool>();
    }
    if (vm.find("hess") != vm.end()) {
        get_hess = vm["hess"].as<bool>();
    }

    // Solve problem:

    try {
        std::ifstream from;
        srs::fopen(from, input_file);

        Gauss_filetype filetype = get_filetype(input_file.c_str());

        Gauss_data gauss(from, filetype);

        natoms = gauss.get_natoms();

        gauss.get_irc_data(mep);
        gauss.get_irc_geom(geom);
        gauss.get_irc_grad(grad);

        if (get_hess) {
            gauss.get_irc_hess(hess);
        }

        if (reverse_mep) {
            srs::dvector mep_tmp(mep.size());
            const int npoints = mep.size() / 2;
            int indx1, indx2;
            for (int i = 0; i < npoints; i++) {
                indx1 = (npoints - 1 - i) * 2;
                indx2 = i * 2;
                for (int j = 0; j < 2; j++) {
                    mep_tmp[indx1 + j] = mep[indx2 + j];
                }
            }
            mep = mep_tmp;
        }
        create_fu31(filetype);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

//-----------------------------------------------------------------------------

Gauss_filetype get_filetype(const char* filename)
{
    std::string suffix = srs::get_suffix(filename);
    if ((suffix == ".out") || (suffix == ".log")) {
        return out;
    }
    else if ((suffix == ".fch") || (suffix == ".fchk")) {
        return fchk;
    }
    else {
        throw Error("unknown suffix: " + suffix);
    }
}

void create_fu31(const Gauss_filetype filetype)
{
    const char* output_file = "gauss2poly.fu31";

    std::ofstream to;
    srs::fopen(to, output_file);

    srs::Format<double> sci;
    srs::Format<double> fix;
    sci.scientific_E().width(16).precision(8);
    fix.fixed().width(12).precision(8);

    double conv_factor;
    if (filetype == out) {
        conv_factor = 1.0;
    }
    else {                         // filetype == fchk
        conv_factor = datum::a_0;  // bohr to angstrom
    }

    const int npoints = mep.size() / 2;
    const int natoms3 = natoms * 3;
    const int nhess   = natoms3 * (natoms3 + 1) / 2;

    double smep;
    double vmep;

    int indx1 = 0;
    int indx2 = 1;
    int indx3 = 0;
    int indx4 = 0;
    int count = 0;

    for (int i = 0; i < npoints; i++) {
        vmep = mep[indx1];
        smep = mep[indx2];
        if (change_sign) {
            smep *= -1.0;
        }
        vmep -= energy0;
        smep += s_corr;

        to << "*POINT \n\n"
           << " SMEP\t" << sci(smep) << "\n\n"
           << " VMEP\t" << sci(vmep) << "\n\n";

        to << " GEOM \n";
        count = 0;
        for (int j = 0; j < natoms3; j++) {
            to << fix(geom[indx3 + j] * conv_factor);
            count++;
            if (count == 3) {
                to << '\n';
                count = 0;
            }
        }
        to << " END \n\n";

        to << " GRADS \n";
        count = 0;
        for (int j = 0; j < natoms3; j++) {
            to << sci(grad[indx3 + j]);
            count++;
            if (count == 3) {
                to << '\n';
                count = 0;
            }
        }
        to << " END \n\n";

        if (get_hess) {
            to << " HESSIANS \n";
            count = 0;
            for (int j = 0; j < nhess; j++) {
                to << sci(hess[indx4 + j]);
                count++;
                if (count == 5) {
                    to << '\n';
                    count = 0;
                }
            }
            if (count != 0) {
                to << '\n';
            }
            to << " END \n\n";
        }
        indx1 += 2;
        indx2 += 2;
        indx3 += natoms3;
        indx4 += nhess;
    }
    std::cout << "output is written to " << output_file << '\n';
}
