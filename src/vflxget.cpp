// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

struct Error : std::runtime_error {
    Error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

double val;
bool get_pressure;
std::string pattern;

//------------------------------------------------------------------------------

void set_values(const std::vector<std::string>& args);
void get_data(std::ifstream& from);

//------------------------------------------------------------------------------

// Program for extracting rate data from VARIFLEX output file at a given
// temperature or pressure.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() < 4) {
        std::cerr << "usage: " << args[0] << " variflex.out value unit [uni]\n";
        return 1;
    }
    std::ifstream from(args[1]);
    if (!from) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }
    try {
        set_values(args);
        get_data(from);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

void set_values(const std::vector<std::string>& args)
{
    val = std::stod(args[2]);
    if (val <= 0.0) {
        throw Error("bad value for T/P: " + std::to_string(val));
    }

    bool get_uni = false;
    pattern = "Pressure  Temp   k_bi-TST";
    if ((args[4] != "") && (args[4] == "uni")) {
        pattern = "Pressure  Temp   k_uni-TST";
        get_uni = true;
    }

    std::string unit = args[3];
    get_pressure = false;
    if ((unit == "Torr") || (unit == "torr")) {
        std::cout << "P = " << val << " " << unit << '\n';
        if (get_uni) {
            std::cout << "T/K\tk_uni-TST/s-1\tk_cid(LowP)/s-1\n";
        }
        else {
            std::cout << "T/K\tk_bi-TST/(cm3/s)\tk_ca(LowP)/(cm3/s)\n";
        }
        get_pressure = true;
    }
    else if ((unit == "Kelvin") || (unit == "kelvin") || (unit == "K")) {
        std::cout << "T = " << val << " " << unit << '\n';
        if (get_uni) {
            std::cout << "P/Torr\tk_uni-TST/s-1\tk_cid(LowP)/s-1\n";
        }
        else {
            std::cout << "P/Torr\tk_bi-TST/(cm3/s)\tk_ca(LowP)/(cm3/s)\n";
        }
    }
    else {
        throw Error("unknown unit: " + unit);
    }
}

void get_data(std::ifstream& from)
{
    std::string line;
    double p, t, k, k0;
    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            std::getline(from, line); // ignore one line
            while (std::getline(from, line)) {
                std::istringstream iss(line);
                iss >> p >> t >> k >> k0;
                if (!iss) {
                    break;
                }
                if (get_pressure) { // get data for given pressure
                    if (p == val) {
                        std::cout
                            << t << '\t'
                            << std::setiosflags(std::ios_base::scientific) << k
                            << "\t\t" << k0 << '\n'
                            << std::resetiosflags(std::ios_base::scientific);
                    }
                }
                else { // get data for given temperature
                    if (t == val) {
                        std::cout
                            << p << '\t'
                            << std::setiosflags(std::ios_base::scientific) << k
                            << "\t\t" << k0 << '\n'
                            << std::resetiosflags(std::ios_base::scientific);
                    }
                }
            }
        }
    }
}

