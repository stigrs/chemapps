////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013-2018 Stig Rune Sellevag. All rights reserved.
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

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4505)  // caused by boost/program_options.hpp>
#endif                           // _MSC_VER

#include <chem/gauss_data.h>
#include <srs/utils.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER


//
//  Summarize NMR shieldings from a Gaussian calculation.
//  The program produces the same output as GaussView.
//
int main(int argc, char* argv[])
{
    // Parse command line arguments:

    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file")
        ("method,m", po::value<std::string>(), "NMR method (default is GIAO)")
        ("tol,t", po::value<double>(), "degeneracy tolerance (default is 0.05)");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    std::string input_file;
    std::string nmr_method = "SCF GIAO";  // default NMR method
    double degen_tol       = 0.05;        // default degeneracy tolerance

    if (vm.find("help") != vm.end()) {
        std::cout << options << '\n';
        return 0;
    }
    if (vm.find("file") != vm.end()) {
        input_file = vm["file"].as<std::string>();
    }
    else {
        std::cerr << options << '\n';
        return 1;
    }
    if (vm.find("method") != vm.end()) {
        nmr_method = vm["method"].as<std::string>();
    }
    if (vm.find("tol") != vm.end()) {
        degen_tol = vm["tol"].as<double>();
    }

    std::ifstream from;
    srs::fopen(from, input_file);

    // Get NMR data:

    std::vector<Gauss_NMR> nmr;

    Gauss_data gauss(from, out);
    gauss.get_nmr_data(nmr, nmr_method, degen_tol);

    // Output results:

    srs::Format<char> hline;
    hline.width(37).fill('-');

    std::cout << "\nSummary of NMR spectrum (" << nmr_method
              << " magnetic shieldings)\n"
              << "Degenerate peaks are condensed together "
              << "(degeneracy tolerance " << degen_tol << ")\n\n"
              << "Shielding/ppm\t"
              << "Degen.\t"
              << "Elem.\t"
              << "Atoms\n"
              << hline('-') << '\n';

    srs::Format<double> fix8;
    fix8.fixed().width(9).precision(4);

    for (auto& ni : nmr) {
        double shield = std::accumulate(ni.shield.begin(), ni.shield.end(), 0.0)
                        / static_cast<double>(ni.shield.size());

        std::cout << fix8(shield) << '\t' << ni.shield.size() << '\t' << ni.atom
                  << '\t';

        std::sort(ni.number.begin(), ni.number.end());

        for (srs::size_t j = 0; j < ni.number.size(); ++j) {
            std::cout << ni.number[j];
            if ((ni.number.size() > 1) && (j < ni.number.size() - 1)) {
                std::cout << ',';
            }
        }
        std::cout << '\n';
    }
}
