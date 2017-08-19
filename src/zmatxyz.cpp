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

#include <chem/molecule.h>
#include <chem/molecule_io.h>
#include <chem/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

//
// Program for converting between XYZ and Z matrix file formats.
//
int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file")
        ("xyz", po::bool_switch()->default_value(false), "convert to XYZ")
        ("zmat", po::bool_switch()->default_value(false), 
         "convert to Z matrix");
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
        chem::fopen(from, input_file);

        Molecule mol(from, to, "Molecule");

        std::string output_file;
        output_file = chem::strip_suffix(input_file, ".inp");

        if (vm["xyz"].as<bool>()) {
            mol.get_zmat()->load(from);
            output_file = output_file + ".xyz";
            chem::fopen(to, output_file.c_str());
            chem::print_xyz_format(to, mol.get_atoms(), mol.get_xyz(), "");
        }
        else if (vm["zmat"].as<bool>()) {
            output_file = output_file + ".zmat";
            chem::fopen(to, output_file.c_str());
            mol.get_zmat()->print(to);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
