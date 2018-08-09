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
#pragma warning(disable : 4458)  // caused by boost/iostreams/filter/gzip.hpp
#pragma warning(disable : 4706)  // caused by boost/iostreams/filter/zlib.hpp
#endif                           // _MSC_VER

#include <chem/gauss_data.h>
#include <chem/ptable.h>
#include <srs/datum.h>
#include <srs/utils.h>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <exception>
#include <fstream>
#include <gsl/gsl>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif  // _MSC_VER

//
// Extracts current Cartesian coordinates from a Gaussian fchk.gz file
// and writes them to standard output in XYZ format.
//
int main(int argc, char* argv[])
{
    auto args = gsl::multi_span<char*>(argv, argc);
    if (argc != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.fchk.gz\n";
        return 1;
    }

    std::string filename       = args[1];
    std::string::size_type pos = filename.rfind(".fchk.gz");
    if (pos == std::string::npos) {
        std::cerr << filename << " is not a fchk.gz file\n";
        return 1;
    }

    // Decompress the fchk.gz file and read the content by reference into
    // a regular istream:

    try {
        std::ifstream file(filename.c_str(),
                           std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(file);
        std::istream from(&in);

        // Get current Cartesian coordinates:

        Gauss_data gauss(from, fchk);
        Gauss_coord coord;

        gauss.get_opt_cart_coord(coord);

        // Write Cartesian coordinates to standard output in XYZ format:

        srs::Format<double> fix8;
        fix8.fixed().width(15).precision(8);

        std::cout << coord.natoms << "\n\n";

        for (int i = 0; i < coord.natoms; ++i) {
            std::cout << ptable::get_atomic_symbol(coord.atnum[i]) << " ";
            for (int j = 0; j < 3; ++j) {
                std::cout << fix8(coord.xyz(i, j) * datum::a_0) << " ";
            }
            std::cout << '\n';
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
    }
}
