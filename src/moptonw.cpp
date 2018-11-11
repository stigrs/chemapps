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
#include <gsl/gsl>
#include <iostream>
#include <stdexcept>
#include <string>


//------------------------------------------------------------------------------

// Error reporting:

struct IO_error : std::runtime_error {
    IO_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

void generate_nw(const std::string& tml_file, const std::string& arc_file);
void get_mopac_geom(const std::string& arc_file);

//------------------------------------------------------------------------------

//
// Program for generating NWChem input file from a template input file and by
// extracting Cartesian coordinates from a MOPAC calculation.
//
int main(int argc, char* argv[])
{
    auto args = gsl::multi_span<char*>(argv, argc);
    if (argc != 3) {
        std::cerr << "Usage: " << args[0] << " nwchem.tml mopac.arc\n\n"
                  << "nwchem.tml: Template file for NWChem input file\n"
                  << "mopac.arc:  Summary file from MOPAC calculation\n";
        return 1;
    }

    try {
        generate_nw(args[1], args[2]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Generate NWChem input file.
void generate_nw(const std::string& tml_file, const std::string& arc_file)
{
    std::ifstream from(tml_file.c_str());
    if (!from) {
        throw IO_error("cannot open " + tml_file);
    }

    const std::string geom_here = "GEOMETRY_HERE";

    std::string line;
    bool found = false;

    while (std::getline(from, line)) {
        if (line.find(geom_here, 0) != std::string::npos) {
            get_mopac_geom(arc_file);
            found = true;
        }
        else {
            std::cout << line << '\n';
        }
    }
    if (!found) {
        throw IO_error("could not find keyword " + geom_here);
    }
}

// Extract optimized Cartesian geometry from MOPAC arc file.
void get_mopac_geom(const std::string& arc_file)
{
    std::ifstream from(arc_file.c_str());
    if (!from) {
        throw IO_error("cannot open " + arc_file);
    }

    const std::string pattern = "FINAL GEOMETRY OBTAINED";

    std::string line;
    std::string atom;
    std::string flag;
    double x;
    double y;
    double z;
    double charge;
    bool found = false;

    srs::Format<double> fix8;
    fix8.fixed().width(15).precision(8);

    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            found = true;
            std::getline(from, line);  // ignore three lines
            std::getline(from, line);
            std::getline(from, line);
            while (from >> atom >> x >> flag >> y >> flag >> z >> flag
                   >> charge) {
                std::cout << "   " << atom << " " << fix8(x) << " " << fix8(y)
                          << " " << fix8(z) << '\n';
            }
        }
    }
    if (!found) {
        throw IO_error("could not find keyword " + pattern);
    }
}
