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

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <chem/gaussian.h>
#include <chem/mopac.h>
#include <chem/gamss.h>
#include <chem/molecule.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program providing Genetic Algorithm Molecular Structure Search.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Genetic Algorithm Molecular Structure Search");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>()) 
		("p,pot", "potential (Gaussian or Mopac)", cxxopts::value<std::string>());
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;
    std::string pot = "Mopac";

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("pot")) {
        pot = args["pot"].as<std::string>();
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

        if (pot == "Gaussian" || pot == "gaussian") {
            Chem::Gamss<Chem::Gaussian> ga(from);
            ga.solve();
        }
        else {
            Chem::Gamss<Chem::Mopac> ga(from);
            ga.solve();
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

