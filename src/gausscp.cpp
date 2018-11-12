////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2018 Stig Rune Sellevag. All rights reserved.
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
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

// Error reporting:

struct IO_error : std::runtime_error {
    IO_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

std::string get_theory(const std::string& filename);
void bsse_se(const std::string& theory, const std::string& filename);

//------------------------------------------------------------------------------

// Program for calculating the BSSE corrected stabilization energy from a
// G03/G09 counterpoise calculation.
//
// Note: Theoretical methods are implemented as needed.
// Current methods: MP2
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "Usage: " << args[0] << " gaussian_file.out\n";
        return 1;
    }

    try {
        std::string theory = get_theory(args[1]);
        bsse_se(theory, args[1]);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

/// Get theoretical method.
std::string get_theory(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    std::string line;
    std::string word;
    std::string theory = "unknown";

    while (std::getline(from, line)) {
        std::istringstream iss(line);
        while (iss >> word) {
            if (std::strcmp(word.substr(0, 1).c_str(), "#") == 0) {
                while (iss >> word) {
                    std::size_t pos = word.find('/');
                    if (pos != std::string::npos) {
                        theory = word.substr(0, pos);
                        break;
                    }
                }
            }
        }
    }
    if (theory == "MP2" || theory == "UMP2") {
        return "EUMP2";
    }
    else {
        throw IO_error(theory + " theory method");
    }
}

// Function for extracting counterpoise BSSE data.
void bsse_se(const std::string& theory, const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    const std::string s_cp_mcbs =
        "Counterpoise: doing MCBS calculation for fragment";
    const std::string s_cp_en = "Counterpoise: corrected energy";
    const std::string s_cp_bsse = "Counterpoise: BSSE energy";

    std::string line;
    std::string ignore;
    std::string word;

    double e_cp = 0.0;    // BSSE corrected energy
    double e_bsse = 0.0;  // BSSE energy
    double e1_mcbs = 0.0; // Energy of fragment 1 using monomer centered basis
    double e2_mcbs = 0.0; // Energy of fragment 2 using momomer centered basis

    while (std::getline(from, line)) {
        std::istringstream iss;
        if (line.find(s_cp_mcbs) != std::string::npos) {
            int frag;
            iss.str(line.substr(s_cp_mcbs.size() + 1));
            iss >> frag;
            while (std::getline(from, line)) {
                std::size_t pos = line.find(theory);
                if (pos != std::string::npos) {
                    iss.clear();
                    iss.str(line.substr(pos + theory.size()));
                    if (frag == 1) {
                        iss >> ignore >> word;
                        e1_mcbs = Stdutils::from_fortran_sci_fmt(word);
                        break;
                    }
                    else if (frag == 2) {
                        iss >> ignore >> word;
                        e2_mcbs = Stdutils::from_fortran_sci_fmt(word);
                        break;
                    }
                    else {
                        throw IO_error("unknown fragment");
                    }
                }
            }
        }
        else if (line.find(s_cp_en) != std::string::npos) {
            iss.clear();
            iss.str(line.substr(s_cp_en.size() + 1));
            iss >> ignore >> e_cp;
        }
        else if (line.find(s_cp_bsse) != std::string::npos) {
            iss.clear();
            iss.str(line.substr(s_cp_bsse.size() + 1));
            iss >> ignore >> e_bsse;
        }
    }
    Stdutils::Format<double> fix8(8);
    Stdutils::Format<double> fix2(2);
    fix8.fixed();
    fix2.fixed();
    std::cout << "\nResults from counterpoise calculation:\n"
              << "--------------------------------------\n"
              << "MCBS energy fragment 1:\t" << fix8(e1_mcbs) << " hartree\n"
              << "MCBS energy fragment 2:\t" << fix8(e2_mcbs) << " hartree\n"
              << "BSSE energy:\t\t" << fix8(e_bsse) << " hartree\n"
              << "BSSE corrected energy:\t" << fix8(e_cp) << " hartree\n"
              << "BSSE corrected stablilization energy:\t"
              << fix2((e_cp - e1_mcbs - e2_mcbs) * 2625.5) << " kJ/mol\n\n"
              << "Data read from " << filename << '\n';
}

