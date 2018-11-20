////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009-2018 Stig Rune Sellevag. All rights reserved.
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
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------

void print_array(std::ostream& to, const std::vector<double>& array);

//-----------------------------------------------------------------------------

const std::string pattern_irc_data =
    "IRC point       1 Results for each geome   R   N=";

const std::string pattern_irc_geom =
    "IRC point       1 Geometries               R   N=";

const std::string pattern_irc_grad =
    "IRC point       1 Gradient at each geome   R   N=";

//-----------------------------------------------------------------------------

// Merge a set of files with Gaussian 98/03 IRC data which have been sorted.
//
// Note: The input files must be in correct order.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() < 2) {
        std::cerr << "usage: " << args[0] << " file1 file2 ... fileN\n";
        return 1;
    }

    try {
        std::vector<double> mep;
        std::vector<double> geom;
        std::vector<double> grad;

        std::ifstream from;
        std::ofstream to;

        const char* output_file = "mergeirc.out";
        Stdutils::fopen(to, output_file);

        int count = 1;
        while (count < argc) {
            auto input_file = args[count];
            std::cout << "Reading " << input_file << " ...\n";
            Stdutils::fopen(from, input_file);

            Chem::Gauss_data gauss(from, Chem::fchk);

            gauss.get_irc_data(mep);
            gauss.get_irc_geom(geom);
            gauss.get_irc_grad(grad);

            from.close();
            count++;
        }

        Stdutils::Format<std::size_t> fmt;
        fmt.width(12);

        to << pattern_irc_data << fmt(mep.size()) << '\n';
        print_array(to, mep);

        to << pattern_irc_geom << fmt(geom.size()) << '\n';
        print_array(to, geom);

        to << pattern_irc_grad << fmt(grad.size()) << '\n';
        print_array(to, grad);

        std::cout << "\nOutput is written to " << output_file << '\n';
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

void print_array(std::ostream& to, const std::vector<double>& array)
{
    Stdutils::Format<double> sci;
    sci.scientific_E().width(16).precision(8);

    int count = 0;
    bool wrote_endl = false;

    for (const auto& vi : array) {
        to << sci(vi);
        count++;
        wrote_endl = false;
        if (count == 5) {
            to << '\n';
            count = 0;
            wrote_endl = true;
        }
    }
    if (!wrote_endl) {
        to << '\n';
    }
}

