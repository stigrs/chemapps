// Copyright (c) 2011-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------

struct Bad_file : std::runtime_error {
    Bad_file(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

const std::string pattern_irc_data =
    "IRC point       1 Results for each geome   R   N=";

const std::string pattern_irc_geom =
    "IRC point       1 Geometries               R   N=";

const std::string pattern_irc_grad =
    "IRC point       1 Gradient at each geome   R   N=";

//------------------------------------------------------------------------------

void print_array(std::ostream& to, const std::vector<double>& array);

//------------------------------------------------------------------------------

// Sort output data from a Gaussian 98/03 IRC calculation for input to
// GaussView so that it displays the IRC walk from reactants to products.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Sort output from Gaussian IRC calculation");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>())
        ("s,sign", "change sign of MEP values (false)", cxxopts::value<bool>()->default_value("false"))
        ("m,rmep", "reverse MEP values (false)", cxxopts::value<bool>()->default_value("false"))
        ("g,rgeom", "reverse geometries and gradients (false)", cxxopts::value<bool>()->default_value("false"))
        ("c,corr", "correction to SMEP (0.0)", cxxopts::value<double>());
    // clang-format on

    auto args = options.parse(argc, argv);

    bool change_sign = false;
    bool reverse_mep = false;
    bool reverse_geom = false;
    double smep_corr = 0.0;
    std::string input_file;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }
    if (args.count("corr")) {
        smep_corr = args["corr"].as<double>();
    }
    if (args.count("sign")) {
        change_sign = args["sign"].as<bool>();
    }
    if (args.count("rmep")) {
        reverse_mep = args["rmep"].as<bool>();
    }
    if (args.count("rgeom")) {
        reverse_geom = args["rgeom"].as<bool>();
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

        Chem::Gauss_filetype filetype;
        std::string suffix = Stdutils::get_suffix(input_file);
        if ((suffix == ".fch") || (suffix == ".fchk")) {
            filetype = Chem::fchk;
        }
        else {
            throw Bad_file("input file is not a fchk file");
        }
        std::string output_file = Stdutils::strip_suffix(input_file, suffix);
        output_file += ".dat";

        Stdutils::fopen(from, input_file);
        Stdutils::fopen(to, output_file);

        // Correct SMEP data and write to output:

        Chem::Gauss_data gauss(from, filetype);

        std::vector<double> mep;
        gauss.get_irc_data(mep);

        if (change_sign) {
            for (std::size_t i = 1; i < mep.size(); i += 2) {
                mep[i] = -1.0 * mep[i] + smep_corr;
            }
        }
        else {
            for (std::size_t i = 1; i < mep.size(); i += 2) {
                mep[i] += smep_corr;
            }
        }

        Stdutils::Format<std::size_t> fmt;
        fmt.width(12);

        to << pattern_irc_data << fmt(mep.size()) << '\n';

        if (reverse_mep) {
            std::vector<double> mep_rev(mep.size());
            const int npoints = static_cast<int>(mep.size() / 2);
            int indx1;
            int indx2;
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

        std::vector<double> geom;
        gauss.get_irc_geom(geom);

        std::vector<double> grad;
        gauss.get_irc_grad(grad);

        if (reverse_geom) {
            const int natoms3 = 3 * gauss.get_natoms();
            const int npoints = static_cast<int>(mep.size() / 2);

            std::vector<double> geom_rev(geom.size());
            std::vector<double> grad_rev(grad.size());

            int indx1;
            int indx2;
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

void print_array(std::ostream& to, const std::vector<double>& array)
{
    Stdutils::Format<double> sci;
    sci.scientific_E().width(16).precision(8);

    int count = 0;
    bool wrote_endl = false;

    for (auto v : array) {
        to << sci(v);
        count++;
        wrote_endl = false;
        if (count == 5) {
            to << '\n';
            count = 0;
            wrote_endl = true;
        }
    }
    if (!wrote_endl) {
        to << '\n';
    }
}

