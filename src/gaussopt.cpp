///////////////////////////////////////////////////////////////////////////////
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
#include <stdutils/stdutils.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>

//  Extracts optimized energies and geometry from Gaussian output file.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.log\n";
        return 1;
    }

    std::ifstream from(args[1]);
    if (!from) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }

    Chem::Gauss_data gauss(from, Chem::out);
    auto en = gauss.get_scf_zpe_energy();
    double en_sum = std::accumulate(en.begin(), en.end(), 0.0);

    Chem::Gauss_coord coord;
    gauss.get_opt_cart_coord(coord);

    Stdutils::Format<double> fix;
    fix.fixed().width(15).precision(8);
    std::cout << "SCF: " << fix(en[0]) << " Hartree\n"
              << "ZPE: " << fix(en[1]) << " Hartree\n"
              << "Tot: " << fix(en_sum) << " Hartree\n\n";

    gauss.print_opt_geom();
}
