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

#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

// Forward declarations:

void convert(const std::string& input_file);

//------------------------------------------------------------------------------

// Converts a XYZ file to GAMESS $DATA format.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " file.xyz\n";
        return 1;
    }

    try {
        convert(args[1]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

//------------------------------------------------------------------------------

void convert(const std::string& input_file)
{
    std::ifstream from;
    Stdutils::fopen(from, input_file);

    std::string line;
    std::string atsymb;
    double x;
    double y;
    double z;

    namespace Pt = Chem::Periodic_table;

    while (from >> atsymb >> x >> y >> z) {
        Stdutils::Format<double> fix(1);
        fix.fixed();
        double atnum = Pt::get_atomic_number(atsymb);
        std::cout << atsymb << "  " << fix(atnum) << "  " << x << "  " << y
                  << "  " << z << '\n';
    }
}

