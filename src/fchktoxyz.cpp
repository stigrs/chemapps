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

#include <chem/gauss_data.h>
#include <chem/periodic_table.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Extracts current Cartesian coordinates from a Gaussian fchk.gz file
// and writes them to standard output in XYZ format.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.fchk.gz\n";
        return 1;
    }

    std::string filename = args[1];
    std::string::size_type pos = filename.rfind(".fchk");
    if (pos == std::string::npos) {
        std::cerr << filename << " is not a fchk file\n";
        return 1;
    }

    try {
        std::ifstream from;
        Stdutils::fopen(from, filename);

        // Get current Cartesian coordinates:

        Chem::Gauss_data gauss(from, Chem::fchk);
        Chem::Gauss_coord coord;

        gauss.get_opt_cart_coord(coord);

        // Write Cartesian coordinates to standard output in XYZ format:

        Stdutils::Format<double> fix8;
        fix8.fixed().width(15).precision(8);

        std::cout << coord.natoms << "\n\n";

        namespace Ptable = Chem::Periodic_table;
        namespace Pc = Numlib::Constants;

        for (int i = 0; i < coord.natoms; ++i) {
            std::cout << Ptable::get_atomic_symbol(coord.atnum[i]) << " ";
            for (int j = 0; j < 3; ++j) {
                std::cout << fix8(coord.xyz(i, j) * Pc::a_0) << " ";
            }
            std::cout << '\n';
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
    }
}

