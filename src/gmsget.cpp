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

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>


//------------------------------------------------------------------------------

// Error reporting:

struct IO_error : std::runtime_error {
    IO_error(std::string s) : std::runtime_error(s) {}
};

struct Run_error : std::runtime_error {
    Run_error(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

void parse_inp(const std::string& inpfile, std::string& runtyp);
void parse_log(const std::string& logfile, const std::string& runtyp);
void extract_optimized_energy(std::istream& from);
void extract_zpe(std::istream& from);

//------------------------------------------------------------------------------

//
// Program for extracting data from GAMESS calculations.
//
// Note: All GAMESS keywords must be in uppercase.
//
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " file.inp file.log\n";
        return 1;
    }

    try {
        std::string runtyp;
        parse_inp(argv[1], runtyp);
        parse_log(argv[2], runtyp);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Get RUNTYP from GAMESS input file (file.inp)
void parse_inp(const std::string& inpfile, std::string& runtyp)
{
    std::ifstream from(inpfile.c_str());
    if (!from) {
        throw IO_error("cannot open " + inpfile);
    }

    std::string line;
    std::string token;

    bool found = false;

    while (std::getline(from, line)) {
        if (line.find("$CONTRL", 0) != std::string::npos) {
            std::istringstream iss(line);
            while (iss >> token) {
                if (token.find("RUNTYP=", 0) != std::string::npos) {
                    runtyp = token.substr(7);
                    found  = true;
                }
            }
        }
    }
    if (!found) {
        throw IO_error("could not find RUNTYP");
    }
}

// Extract data from GAMESS output file.
void parse_log(const std::string& logfile, const std::string& runtyp)
{
    std::ifstream from(logfile.c_str());
    if (!from) {
        throw IO_error("cannot open " + logfile);
    }
    if (runtyp == "OPTIMIZE") {
        extract_optimized_energy(from);
    }
    else if (runtyp == "HESSIAN") {
        extract_zpe(from);
    }
    else {
        throw Run_error("unknown RUNTYP=" + runtyp);
    }
}

// Extract electronic energy for optimized geometry.
void extract_optimized_energy(std::istream& from)
{
    const std::string pat_opt = "***** EQUILIBRIUM GEOMETRY LOCATED *****";
    const std::string pat_e   = "TOTAL ENERGY";

    std::string line;
    std::string tmp;
    std::string energy;  // a bit dangerous, but easy (avoids loss of decimals)

    while (std::getline(from, line)) {
        if (line.find(pat_opt, 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find(pat_e, 0) != std::string::npos) {
                    std::istringstream iss(line);
                    iss >> tmp >> tmp >> tmp >> energy;
                    std::cout << "  " << pat_e << " = " << energy << '\n';
                    return;  // success
                }
            }
        }
    }
    throw IO_error("extracting optimized energy failed");
}

// Extract zero-point vibrational energy.
void extract_zpe(std::istream& from)
{
    const std::string pat_zpe = "THE HARMONIC ZERO POINT ENERGY IS";

    std::string line;
    std::string zpe;  // a bit dangerous, but easy (avoids loss of decimals)

    while (std::getline(from, line)) {
        if (line.find(pat_zpe, 0) != std::string::npos) {
            std::getline(from, line);
            std::istringstream iss(line);
            iss >> zpe;
            std::cout << "  ZERO POINT ENERGY = " << zpe << '\n';
            return;  // success
        }
    }
    throw IO_error("extracting zero-point energy failed");
}
