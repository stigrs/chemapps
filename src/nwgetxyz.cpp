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

#include <stdutils/stdutils.h>
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

void extract_geometry(const std::string& outfile);

//------------------------------------------------------------------------------

// Program for extracting optimized geometry from NWChem calculations.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "Usage: " << args[0] << " file.out\n";
        return 1;
    }

    try {
        extract_geometry(args[1]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Extract geometry from NWChem output file.
void extract_geometry(const std::string& outfile)
{
    std::ifstream from(outfile.c_str());
    if (!from) {
        throw IO_error("cannot open " + outfile);
    }

    const std::string pat_opt = "Optimization converged";
    const std::string pat_geom = "No.       Tag          Charge";

    std::string line;
    std::string word;
    int i;
    double charge;
    double x;
    double y;
    double z;

    bool found = false;

    Stdutils::Format<double> fix1;
    Stdutils::Format<double> fix8;
    fix1.fixed().width(5).precision(1);
    fix8.fixed().width(15).precision(8);

    while (std::getline(from, line)) {
        if (line.find(pat_opt, 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find(pat_geom, 0) != std::string::npos) {
                    found = true;
                    std::getline(from, line); // ignore one line
                    while (from >> i >> word >> charge >> x >> y >> z) {
                        std::cout << word << '\t' << fix1(charge) << " "
                                  << fix8(x) << " " << fix8(y) << " " << fix8(z)
                                  << '\n';
                    }
                }
            }
        }
    }
    if (!found) {
        throw IO_error("could not find optimized geometry");
    }
}

