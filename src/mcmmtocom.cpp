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
#pragma warning(disable : 4100 4505)  // caused by boost/program_options.hpp
#endif                                // _MSC_VER

#include <srs/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER


//
// Program for generating Gaussian input files from MCMM solver output.
//
int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file") 
		("proc,N", po::value<int>(), "number of processors") 
		("charge,c", po::value<int>(), "charge of molecule") 
		("spin,s", po::value<int>(), "spin multiplicity of molecule") 
        ("key,k", po::value<std::string>(), "Gaussian keywords") 
        ("title,t", po::value<std::string>(), "title line");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    std::string input_file;
    std::string keywords = "opt freq hf/sto-3g";
    std::string title    = "Title";
    int nproc            = 1;
    int charge           = 0;
    int spin             = 1;

    if (vm.find("help") != vm.end()) {
        std::cout << options << '\n';
        return 0;
    }
    if (vm.find("title") != vm.end()) {
        title = vm["title"].as<std::string>();
    }
    if (vm.find("key") != vm.end()) {
        keywords = vm["key"].as<std::string>();
    }
    if (vm.find("proc") != vm.end()) {
        nproc = vm["proc"].as<int>();
    }
    if (vm.find("charge") != vm.end()) {
        charge = vm["charge"].as<int>();
    }
    if (vm.find("spin") != vm.end()) {
        spin = vm["spin"].as<int>();
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
        srs::fopen(from, input_file);
        std::string line;
        while (std::getline(from, line)) {
            if (line.find("Conformer:") != std::string::npos) {
                std::string token;
                int nc;
                std::istringstream iss(line);
                iss >> token >> nc;

                std::ofstream to;
                std::string base = srs::strip_suffix(input_file, ".out");
                base += "_c" + srs::to_string(nc);
                std::string com = base + ".com";
                srs::fopen(to, com.c_str());

                to << "%nprocshared=" << nproc << '\n'
                   << "%chk=" << base << ".chk\n"
                   << "# " << keywords << "\n\n"
                   << title << "\n\n"
                   << charge << " " << spin << '\n';

                for (int it = 0; it < 5; ++it) {
                    std::getline(from, line);  // ignore
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
                        srs::Format<double> fix;
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
