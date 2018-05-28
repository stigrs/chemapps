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

#include <srs/utils.h>
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

//------------------------------------------------------------------------------

// Forward declarations:

void extract_geometry(const std::string& logfile);

//------------------------------------------------------------------------------

//
//   Program for extracting optimized geometry from GAMESS calculations.
//
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " file.log\n";
        return 1;
    }

    try {
        extract_geometry(argv[1]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Extract geometry from GAMESS output file.
void extract_geometry(const std::string& logfile)
{
    std::ifstream from(logfile.c_str());
    if (!from) {
        throw IO_error("cannot open " + logfile);
    }

    const std::string pat_opt  = "***** EQUILIBRIUM GEOMETRY LOCATED *****";
    const std::string pat_geom = "COORDINATES OF ALL ATOMS ARE";

    std::string line;
    std::string atom;
    double charge;
    double x;
    double y;
    double z;

    bool found = false;

    srs::Format<double> fix1;
    srs::Format<double> fix10;
    fix1.fixed().width(5).precision(1);
    fix10.fixed().width(15).precision(10);

    while (std::getline(from, line)) {
        if (line.find(pat_opt, 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find(pat_geom, 0) != std::string::npos) {
                    found = true;
                    std::getline(from, line);  // ignore two lines
                    std::getline(from, line);
                    while (from >> atom >> charge >> x >> y >> z) {
                        std::cout << atom << '\t' << fix1(charge) << " "
                                  << fix10(x) << " " << fix10(y) << " "
                                  << fix10(z) << '\n';
                    }
                }
            }
        }
    }
    if (!found) {
        throw IO_error("could not find optimized geometry");
    }
}
