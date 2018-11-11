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
#include <gsl/gsl>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


//------------------------------------------------------------------------------

// Error reporting:

struct IO_error : std::runtime_error {
    IO_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

void parse_nw(const std::string& nwfile,
              std::vector<std::pair<std::string, std::string>>& task);
void parse_out(const std::string& outfile,
               std::vector<std::pair<std::string, std::string>>& task);
void extract_energy(std::istream& from, const std::string& theory);
void extract_optimized_energy(std::istream& from, const std::string& theory);
void extract_zpe(std::istream& from, const std::string& theory);

//------------------------------------------------------------------------------

//
// Program for extracting data from NWChem calculations.
//
int main(int argc, char* argv[])
{
    auto args = gsl::multi_span<char*>(argv, argc);
    if (argc != 3) {
        std::cout << "Usage: " << args[0] << " file.nw file.out\n";
        return 1;
    }

    std::vector<std::pair<std::string, std::string>> task;

    try {
        parse_nw(args[1], task);
        parse_out(args[2], task);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Get TASK directives from NWChem input file (file.nw)
void parse_nw(const std::string& nwfile,
              std::vector<std::pair<std::string, std::string>>& task)
{
    std::ifstream from(nwfile.c_str());
    if (!from) {
        throw IO_error("cannot open " + nwfile);
    }

    std::string line;
    std::string token;
    std::string theory;
    std::string operation;

    std::pair<std::string, std::string> task_tmp;

    bool found = false;

    while (std::getline(from, line)) {
        if (line.find("task", 0) != std::string::npos) {
            found = true;
            std::istringstream iss(line);
            iss >> token >> theory >> operation;
            task_tmp.first  = theory;
            task_tmp.second = operation;
            task.push_back(task_tmp);
        }
    }
    if (!found) {
        throw IO_error("could not find TASK directive");
    }
}

// Extract data from NWChem output file.
void parse_out(const std::string& outfile,
               std::vector<std::pair<std::string, std::string>>& task)
{
    std::ifstream from(outfile.c_str());
    if (!from) {
        throw IO_error("cannot open " + outfile);
    }

    for (auto it : task) {
        std::cout << "TASK " << it.first << ' ' << it.second << ":\n";
        if (it.second == "optimize") {
            extract_optimized_energy(from, it.first);
        }
        else if (it.second == "freq") {
            extract_zpe(from, it.first);
        }
        else {
            extract_energy(from, it.first);
        }
    }
}

// Extract single-point electronic energy.
void extract_energy(std::istream& from, const std::string& theory)
{
    std::string pat_mod;
    std::string pat_e;

    if (theory == "scf") {
        pat_mod = "NWChem SCF Module";
        pat_e   = "Total SCF energy";
    }
    else if (theory == "dft") {
        pat_mod = "NWChem DFT Module";
        pat_e   = "Total DFT energy";
    }

    std::string line;
    std::string tmp;
    std::string energy;  // a bit dangerous, but easy (avoids loss of decimals)

    while (std::getline(from, line)) {
        if (line.find(pat_mod, 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find(pat_e, 0) != std::string::npos) {
                    std::istringstream iss(line);
                    iss >> tmp >> tmp >> tmp >> tmp >> energy;
                    std::cout << "  " << pat_e << " = " << energy << '\n';
                    return;  // success
                }
            }
        }
    }
    throw IO_error("extracting optimized energy failed");
}

// Extract electronic energy for optimized geometry.
void extract_optimized_energy(std::istream& from, const std::string& theory)
{
    std::string pat_opt;
    std::string pat_mod;
    std::string pat_conv;
    std::string pat_e;

    pat_opt  = "NWChem Geometry Optimization";
    pat_conv = "Optimization converged";

    if (theory == "scf") {
        pat_mod = "NWChem SCF Module";
        pat_e   = "Total SCF energy";
    }
    else if (theory == "dft") {
        pat_mod = "NWChem DFT Module";
        pat_e   = "Total DFT energy";
    }

    std::string line;
    std::string tmp;
    std::string energy;  // a bit dangerous, but easy (avoids loss of decimals)

    while (std::getline(from, line)) {
        if (line.find(pat_opt, 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find(pat_mod, 0) != std::string::npos) {
                    while (std::getline(from, line)) {
                        if (line.find(pat_e, 0) != std::string::npos) {
                            std::istringstream iss(line);
                            iss >> tmp >> tmp >> tmp >> tmp >> energy;
                        }
                        if (line.find(pat_conv, 0) != std::string::npos) {
                            std::cout << "  " << pat_e << " = " << energy
                                      << '\n';
                            return;  // success
                        }
                    }
                }
            }
        }
    }
    throw IO_error("extracting optimized energy failed");
}

// Extract zero-point vibrational energy.
void extract_zpe(std::istream& from, const std::string& theory)
{
    std::string pat_freq;
    std::string pat_zpe;
    std::string pat_mod;

    pat_freq = "NWChem Nuclear Hessian and Frequency Analysis";
    pat_zpe  = "Zero-Point correction to Energy";

    if (theory == "scf") {
        pat_mod = "NWChem SCF Module";
    }
    else if (theory == "dft") {
        pat_mod = "NWChem DFT Module";
    }

    std::string line;
    std::string tmp;
    std::string zpe;  // a bit dangerous, but easy (avoids loss of decimals)

    while (std::getline(from, line)) {
        if (line.find(pat_freq, 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find(pat_mod, 0) != std::string::npos) {
                    while (std::getline(from, line)) {
                        if (line.find(pat_zpe, 0) != std::string::npos) {
                            std::istringstream iss(line);
                            iss >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp
                                >> tmp >> zpe;
                            std::cout << "  " << pat_zpe << " = " << zpe
                                      << '\n';
                            return;  // success
                        }
                    }
                }
            }
        }
    }
    throw IO_error("extracting zero-point energy failed");
}
