////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011-2018 Stig Rune Sellevag. All rights reserved.
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

struct Bad_file : std::runtime_error {
    Bad_file(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

const char pattern_irc_data[]
    = "IRC point       1 Results for each geome   R   N=";

const char pattern_irc_geom[]
    = "IRC point       1 Geometries               R   N=";

const char pattern_irc_grad[]
    = "IRC point       1 Gradient at each geome   R   N=";

//------------------------------------------------------------------------------

void print_array(std::ostream& to, const srs::dvector& array);

//------------------------------------------------------------------------------

//
// Sort output data from a Gaussian 98/03 IRC calculation for input to
// GaussView so that it displays the IRC walk from reactants to products.
//
int main(int argc, char* argv[])
{
    // Parse command line options:

    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file")
        ("sign,cs", po::bool_switch()->default_value(false), "change sign of MEP values (false)")
        ("rmep,rs", po::bool_switch()->default_value(false), "reverse MEP values (false)")
        ("rgeom,rx", po::bool_switch()->default_value(false), "reverse geometries and gradients (false)")
        ("corr,c", po::value<double>(), "correction to SMEP (0.0)");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    bool change_sign  = false;
    bool reverse_mep  = false;
    bool reverse_geom = false;
    double smep_corr  = 0.0;
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
    if (vm.find("corr") != vm.end()) {
        smep_corr = vm["corr"].as<double>();
    }
    if (vm.find("sign") != vm.end()) {
        change_sign = vm["sign"].as<bool>();
    }
    if (vm.find("rmep") != vm.end()) {
        reverse_mep = vm["rmep"].as<bool>();
    }
    if (vm.find("rgeom") != vm.end()) {
        reverse_geom = vm["rgeom"].as<bool>();
    }

    if (change_sign) {
        std::cout << "Change sign of SMEP:\tyes\n";
    }
    else {
        std::cout << "Change sign of SMEP:\tno\n";
    }
    if (reverse_mep) {
        std::cout << "Reverse MEP data:\tyes\n";
    }
    else {
        std::cout << "Reverse MEP data:\tno\n";
    }
    if (reverse_geom) {
        std::cout << "Reverse geom./grad.:\tyes\n";
    }
    else {
        std::cout << "Reverse geom./grad.:\tno\n";
    }
    std::cout << "SMEP correction:\t" << smep_corr << '\n';


    try {
        std::ifstream from;
        std::ofstream to;

        Gauss_filetype filetype;
        std::string suffix = srs::get_suffix(input_file);
        if ((suffix == ".fch") || (suffix == ".fchk")) {
            filetype = fchk;
        }
        else {
            throw Bad_file("input file is not a fchk file");
        }
        std::string output_file = srs::strip_suffix(input_file, suffix);
        output_file += ".dat";

        srs::fopen(from, input_file);
        srs::fopen(to, output_file.c_str());

        // Correct SMEP data and write to output:

        Gauss_data gauss(from, filetype);

        srs::dvector mep;
        gauss.get_irc_data(mep);

        if (change_sign) {
            for (srs::size_t i = 1; i < mep.size(); i += 2) {
                mep[i] = -1.0 * mep[i] + smep_corr;
            }
        }
        else {
            for (srs::size_t i = 1; i < mep.size(); i += 2) {
                mep[i] += smep_corr;
            }
        }

        srs::Format<srs::size_t> fmt;
        fmt.width(12);

        to << pattern_irc_data << fmt(mep.size()) << '\n';

        if (reverse_mep) {
            srs::dvector mep_rev(mep.size());
            const int npoints = mep.size() / 2;
            int indx1, indx2;
            for (int i = 0; i < npoints; i++) {
                indx1 = (npoints - 1 - i) * 2;
                indx2 = i * 2;
                for (int j = 0; j < 2; j++) {
                    mep_rev[indx1 + j] = mep[indx2 + j];
                }
            }
            print_array(to, mep_rev);
        }
        else {
            print_array(to, mep);
        }

        // Correct geometries/gradients and write to output:

        srs::dvector geom;
        gauss.get_irc_geom(geom);

        srs::dvector grad;
        gauss.get_irc_grad(grad);

        if (reverse_geom) {
            const int natoms3 = 3 * gauss.get_natoms();
            const int npoints = mep.size() / 2;

            srs::dvector geom_rev(geom.size());
            srs::dvector grad_rev(grad.size());

            int indx1, indx2;
            for (int i = 0; i < npoints; i++) {
                indx1 = (npoints - 1 - i) * natoms3;
                indx2 = i * natoms3;
                for (int j = 0; j < natoms3; j++) {
                    geom_rev[indx1 + j] = geom[indx2 + j];
                    grad_rev[indx1 + j] = grad[indx2 + j];
                }
            }

            to << pattern_irc_geom << fmt(geom_rev.size()) << '\n';
            print_array(to, geom_rev);

            to << pattern_irc_grad << fmt(grad_rev.size()) << '\n';
            print_array(to, grad_rev);
        }
        else {
            to << pattern_irc_geom << fmt(geom.size()) << '\n';
            print_array(to, geom);

            to << pattern_irc_grad << fmt(grad.size()) << '\n';
            print_array(to, grad);
        }
        std::cout << "Output is written to " << output_file << '\n';
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

void print_array(std::ostream& to, const srs::dvector& array)
{
    srs::Format<double> sci;
    sci.scientific_E().width(16).precision(8);

    int count       = 0;
    bool wrote_endl = false;

    for (srs::size_t i = 0; i < array.size(); i++) {
        to << sci(array[i]);
        count++;
        wrote_endl = false;
        if (count == 5) {
            to << '\n';
            count      = 0;
            wrote_endl = true;
        }
    }
    if (!wrote_endl) {
        to << '\n';
    }
}
