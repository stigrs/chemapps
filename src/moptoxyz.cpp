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

#include <srs/utils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


//------------------------------------------------------------------------------

// Error reporting:

struct IO_error : std::runtime_error {
    IO_error(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

void get_xyz(const std::string& arc_file);

//------------------------------------------------------------------------------

//
// Program for extracting XYZ coordinates from a MOPAC calculation.
//
int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " mopac.arc\n\n"
                  << "mopac.arc:  Summary file from MOPAC calculation\n";
        return 1;
    }

    try {
        get_xyz(argv[1]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Extract optimized Cartesian geometry from MOPAC arc file.
void get_xyz(const std::string& arc_file)
{
    std::ifstream from(arc_file.c_str());
    if (!from) {
        throw IO_error("cannot open " + arc_file);
    }

    // Get data:

    const std::string pattern = "FINAL GEOMETRY OBTAINED";

    std::vector<std::string> atoms;
    std::vector<double> xyz;

    std::string title_line;
    std::string line;
    std::string element;
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
            std::getline(from, line);  // ignore this line
            std::getline(from, line);
            title_line = line;
            std::getline(from, line);  // ignore this line
            while (from >> element >> x >> flag >> y >> flag >> z >> flag
                   >> charge) {
                atoms.push_back(element);
                xyz.push_back(x);
                xyz.push_back(y);
                xyz.push_back(z);
            }
        }
    }
    if (!found) {
        throw IO_error("could not find keyword " + pattern);
    }

    // Write data:

    std::cout << atoms.size() << '\n' << title_line << '\n';

    for (unsigned i = 0; i < atoms.size(); ++i) {
        std::cout << atoms[i] << " ";
        for (unsigned j = 0; j < 3; ++j) {
            std::cout << fix8(xyz[i * 3 + j]) << " ";
        }
        std::cout << '\n';
    }
}
