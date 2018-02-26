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

#include <chem/collision.h>
#include <srs/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

//
// Program for computing the biased random walk model of R. G. Gilbert as
// published in the papers:
//
// "Collision energy exchange in highly vibrationally excited molecules:
// The biased random walk model". J. Chem. Phys., 1984, vol. 80, p. 5501.
//
// "Modeling collisional energy transfer in highly excited molecules"
// J. Chem. Phys., 1990, vol. 92, p. 1819.
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

        Collision model(from);
        model.biased_random_walk(to);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
