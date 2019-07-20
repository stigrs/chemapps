// Copyright (c) 2019 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------
//
// Global variables:

bool avg_shield;
bool avg_sscc;

int atoms;
int center;

std::vector<int> nuclei;

//------------------------------------------------------------------------------
//
// Forward declarations:

void gaussian(std::istream& from);
double gauss_sscc(std::istream& from, const std::string& search_str);
double average(const std::vector<double>& values);

//------------------------------------------------------------------------------
//
// Program for averaging NMR absolute shieldings and indirect spin-spin
// coupling constants from quantum chemistry calculations.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Average NMR shieldings and indirect spin-spin coupling constants");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>()) 
        ("N,atoms", "number of atoms", cxxopts::value<int>())
        ("c,center", "atomic center", cxxopts::value<int>())
        ("n,nuclei", "nucleis to average", cxxopts::value<std::string>())
		("p,prog", "quantum chemistry program", cxxopts::value<std::string>())
        ("s,shield", "average shieldings", cxxopts::value<bool>()->default_value("false"))
        ("j,sscc", "average indirect SSCC", cxxopts::value<bool>()->default_value("true"));

    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;
    std::string prog = "Gaussian";

    avg_shield = false;
    avg_sscc = true;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("prog")) {
        prog = args["prog"].as<std::string>();
    }
    if (args.count("shield")) {
        avg_shield = args["shield"].as<bool>();
    }
    if (args.count("sscc")) {
        avg_sscc = args["sscc"].as<bool>();
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }
    if (args.count("atoms")) {
        atoms = args["atoms"].as<int>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }
    if (args.count("center")) {
        center = args["center"].as<int>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }
    if (args.count("nuclei")) {
        auto tmp = args["nuclei"].as<std::string>();
        int i;
        std::istringstream iss(tmp);
        while (iss >> i) {
            nuclei.push_back(i);
        }
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        Stdutils::fopen(from, input_file);

        if (prog == "Gaussian" || prog == "gaussian") {
            gaussian(from);
        }
        else {
            std::cerr << "Only Gaussian is currently supported\n";
            return 1;
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

void gaussian(std::istream& from)
{
    if (avg_sscc) {
        std::string str = "Fermi Contact (FC) contribution to J";
        std::cout << "J_FC =  " << gauss_sscc(from, str) << '\n';

        str = "Spin-dipolar (SD) contribution to J";
        std::cout << "J_SD =  " << gauss_sscc(from, str) << '\n';

        str = "Paramagnetic spin-orbit (PSO) contribution to J";
        std::cout << "J_PSO = " << gauss_sscc(from, str) << '\n';

        str = "Diamagnetic spin-orbit (DSO) contribution to J";
        std::cout << "J_DSO = " << gauss_sscc(from, str) << '\n';

        str = "Total nuclear spin-spin coupling J";
        std::cout << "J_tot = " << gauss_sscc(from, str) << '\n';
    }
}

double gauss_sscc(std::istream& from, const std::string& search_str)
{
    from.clear();
    from.seekg(0, std::ios_base::beg); // move to beginning of file

    std::string line;
    std::string word;

    int ibuf;
    int pos_center = 0;
    int nuc_start = 0;

    bool center_found = false;

    std::vector<double> sscc;

    while (std::getline(from, line)) {
        if (line.find(search_str) != std::string::npos) {
            while (std::getline(from, line)) {
                std::istringstream iss(line);
                while (iss >> ibuf) {
                    ++pos_center;
                    if (ibuf == center) {
                        center_found = true;
                        break;
                    }
                }
                if (center_found) {
                    break;
                }
                else {
                    pos_center = 0;
                    for (int it = nuc_start; it < atoms; ++it) {
                        std::getline(from, line); // ignore
                    }
                    if (atoms > 5) {
                        nuc_start += 5;
                    }
                }
            }
            if (center_found) {
                int inuc = 0;
                for (int it = nuc_start; it < atoms; ++it) {
                    std::getline(from, line);
                    std::istringstream iss(line);
                    iss >> ibuf;
                    if (ibuf == nuclei[inuc]) {
                        ++inuc;
                        for (int i = 0; i < pos_center; ++i) {
                            iss >> word;
                        }
                        sscc.push_back(Stdutils::from_fortran_sci_fmt(word));
                    }
                }
                break;
            }
        }
    }
    return average(sscc);
}

double average(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

