////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

// Error reporting:

struct File_error : std::runtime_error {
    File_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

void convert(const std::string& input_file, const std::string& fmt);

//------------------------------------------------------------------------------

//
// Converts a Materials Studio CAR file to XYZ format.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() < 2) {
        std::cerr << "usage: " << args[0] << " file.car [format]\n\n"
                  << "format = gaussian (default)\n"
                  << "         gamess\n";
        return 1;
    }

    try {
        std::string fmt = "gaussian";
        if (argc == 3) {
            fmt = args[2];
        }
        convert(args[1], fmt);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

//------------------------------------------------------------------------------

void convert(const std::string& input_file, const std::string& fmt)
{
    std::ifstream from;
    Stdutils::fopen(from, input_file);

    std::string line, label, atsymb, s1, s2;
    double x, y, z, dd;
    int ii;

    for (int i = 0; i < 4; ++i) { // ignore four first lines
        if (!std::getline(from, line)) {
            throw File_error(input_file + " is corrupt");
        }
    }

    while (std::getline(from, line)) {
        if (line.find("end") != std::string::npos) {
            break;
        }
        std::istringstream iss(line);
        iss >> label >> x >> y >> z >> s1 >> ii >> s2 >> atsymb >> dd;
        if (!iss) {
            throw File_error(input_file + " is corrupt");
        }
        if (fmt == "gamess") {
            Stdutils::Format<double> fix(1);
            fix.fixed();
            double atnum = Chem::Periodic_table::get_atomic_number(atsymb);
            std::cout << atsymb << "  " << fix(atnum) << "  " << x << "  " << y
                      << "  " << z << '\n';
        }
        else { // Gaussian
            std::cout << atsymb << "  " << x << "  " << y << "  " << z << '\n';
        }
    }
}

