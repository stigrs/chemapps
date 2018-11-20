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

#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Program for generating Gaussian input files from MCMM solver output.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Generate Gaussian input files from MCMM solver output");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>()) 
		("proc,N", "number of processors", cxxopts::value<int>()) 
		("charge,c", "net molecular charge", cxxopts::value<int>()) 
		("spin,s", "spin multiplicity", cxxopts::value<int>()) 
        ("key,k", "Gaussian keywords", cxxopts::value<std::string>()) 
        ("title,t", "title line", cxxopts::value<std::string>());
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;
    std::string keywords = "opt freq hf/sto-3g";
    std::string title = "Title";
    int nproc = 1;
    int charge = 0;
    int spin = 1;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("title")) {
        title = args["title"].as<std::string>();
    }
    if (args.count("key")) {
        keywords = args["key"].as<std::string>();
    }
    if (args.count("proc")) {
        nproc = args["proc"].as<int>();
    }
    if (args.count("charge")) {
        charge = args["charge"].as<int>();
    }
    if (args.count("spin")) {
        spin = args["spin"].as<int>();
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
        std::string line;
        while (std::getline(from, line)) {
            if (line.find("Conformer:") != std::string::npos) {
                std::string token;
                int nc;
                std::istringstream iss(line);
                iss >> token >> nc;

                std::ofstream to;
                std::string base = Stdutils::strip_suffix(input_file, ".out");
                base += "_c" + std::to_string(nc);
                std::string com = base + ".com";
                Stdutils::fopen(to, com.c_str());

                to << "%nprocshared=" << nproc << '\n'
                   << "%chk=" << base << ".chk\n"
                   << "# " << keywords << "\n\n"
                   << title << "\n\n"
                   << charge << " " << spin << '\n';

                for (int it = 0; it < 5; ++it) {
                    std::getline(from, line); // ignore
                }
                int center;
                std::string atom;
                double x;
                double y;
                double z;
                while (std::getline(from, line)) {
                    iss = std::istringstream(line);
                    if (line.find("--") != std::string::npos) {
                        break;
                    }
                    if (iss >> center >> atom >> x >> y >> z) {
                        Stdutils::Format<double> fix;
                        fix.fixed().width(10).precision(6);
                        to << atom << '\t' << fix(x) << "  " << fix(y) << "  "
                           << fix(z) << '\n';
                    }
                }
                to << '\n';
            }
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

