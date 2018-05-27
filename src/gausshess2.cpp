///////////////////////////////////////////////////////////////////////////////
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

#include <srs/utils.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


//
//  Extracts Cartesian force constants from Gaussian output file.
//
int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " gaussian.log\n";
        return 1;
    }
    std::ifstream from(argv[1]);
    if (!from) {
        std::cerr << "cannot open " << argv[1] << '\n';
        return 1;
    }

    const char pattern[] = "Force constants in Cartesian coordinates:";

    bool found = false;
    double fc;
    int token;
    std::string line, data, fc_str;
    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            found     = true;
            int count = 0;
            std::cout << " HESSIAN\n";
            while (from >> token) {
                std::getline(from, data);
                std::istringstream iss(data);
                while (iss >> fc_str) {  // Not fool-proof, be careful!
                    fc = srs::from_fortran_sci_fmt(fc_str);
                    if (static_cast<int>(fc) / fc != 1.0) {
                        std::cout
                            << std::setw(16) << std::setprecision(8)
                            << std::setiosflags(std::ios_base::scientific
                                                | std::ios_base::uppercase)
                            << fc;
                    }
                }
            }
            if (count != 0) {
                std::cout << '\n';
            }
            std::cout << " END\n\n";
        }
    }
    if (!found) {
        std::cerr << "could not find force constants\n";
        return 1;
    }
}
