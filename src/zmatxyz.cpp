////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <chem/molecule.h>
#include <chem/impl/io_support.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

// Program for converting between XYZ and Z matrix file formats.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Convert between XYZ and Z-matrix file formats");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>())
        ("xyz", "convert to XYZ", cxxopts::value<bool>()->default_value("false"))
        ("zmat", "convert to Z-matrix", cxxopts::value<bool>()->default_value("false"));
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        Stdutils::fopen(from, input_file);

        Chem::Molecule mol(from);

        std::ofstream to;
        std::string output_file;
        output_file = Stdutils::strip_suffix(input_file, ".inp");

        if (args["xyz"].as<bool>()) {
            mol.int_coord().load(from);
            output_file = output_file + ".xyz";
            Stdutils::fopen(to, output_file.c_str());
            Chem::Impl::print_xyz_format(to, mol.atoms(), mol.cart_coord(), "");
        }
        else if (args["zmat"].as<bool>()) {
            output_file = output_file + ".zmat";
            Stdutils::fopen(to, output_file.c_str());
            mol.int_coord().print(to);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

