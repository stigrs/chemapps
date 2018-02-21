//////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
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
///////////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4100 4505)  // caused by boost/program_options.hpp
#endif                                // _MSC_VER

#include <chem/mcmm.h>
#include <chem/molecule.h>
#include <chem/mopac.h>
#include <srs/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER


//
// Program providing Monte Carlo Multiple Minima (MCMM) solver.
//
int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    std::string input_file;

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

    try {
        std::ifstream from;
        std::ofstream to;

        std::string output_file;
        output_file = srs::strip_suffix(input_file, ".inp");
        output_file = output_file + ".out";

        srs::fopen(from, input_file);
        srs::fopen(to, output_file.c_str());

        Molecule mol(from, to);
        Mcmm<Mopac> mc(from, mol, "Mcmm", true);
        mc.solve(to);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
