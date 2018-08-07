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

#include <srs/utils.h>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

struct Error : std::runtime_error {
    Error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

double val;
bool get_pressure;
std::string pattern;

//------------------------------------------------------------------------------

void set_values(char* argv[]);
void get_data(std::ifstream& from);

//------------------------------------------------------------------------------

//
// Program for extracting rate data from VARIFLEX output file at a given
// temperature or pressure.
//
int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "usage: " << argv[0] << " variflex.out value unit [uni]\n";
        return 1;
    }
    std::ifstream from(argv[1]);
    if (!from) {
        std::cerr << "cannot open " << argv[1] << '\n';
        return 1;
    }
    try {
        set_values(argv);
        get_data(from);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

void set_values(char* argv[])
{
    val = std::atof(argv[2]);
    if (val <= 0.0) {
        throw Error("bad value for T/P: " + srs::to_string(val));
    }

    bool get_uni = false;
    pattern      = "Pressure  Temp   k_bi-TST";
    if ((argv[4] != 0) && (std::strcmp(argv[4], "uni") == 0)) {
        pattern = "Pressure  Temp   k_uni-TST";
        get_uni = true;
    }

    std::string unit = argv[3];
    get_pressure     = false;
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
            std::getline(from, line);  // ignore one line
            while (std::getline(from, line)) {
                std::istringstream iss(line);
                iss >> p >> t >> k >> k0;
                if (!iss) {
                    break;
                }
                if (get_pressure) {  // get data for given pressure
                    if (p == val) {
                        std::cout
                            << t << '\t'
                            << std::setiosflags(std::ios_base::scientific) << k
                            << "\t\t" << k0 << '\n'
                            << std::resetiosflags(std::ios_base::scientific);
                    }
                }
                else {  // get data for given temperature
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
