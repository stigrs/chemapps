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

#include <chem/gauss_data.h>
#include <srs/array.h>
#include <srs/utils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

//
// Program for extracting data from a Gaussian relaxed PES scan.
//
int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " gaussian_output_file\n";
        return 1;
    }

    try {
        std::ifstream from;
        srs::fopen(from, argv[1]);

        std::string scan_coord;
        srs::dvector coord;
        srs::dvector energy;

        Gauss_data gauss(from, out);
        gauss.get_pes_scan_data(scan_coord, coord, energy);

        srs::Format<double> fix6;
        fix6.fixed();

        std::cout << "#\t" << scan_coord << "\t\t"
                  << "Eigenvalues\n";
        for (int i = 0; i < energy.size(); ++i) {
            std::cout << i + 1 << '\t' << fix6(coord[i]) << '\t'
                      << fix6(energy[i]) << '\n';
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
