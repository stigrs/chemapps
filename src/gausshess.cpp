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

#include <stdutils/stdutils.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

//  Extracts Cartesian force constants from Gaussian fchk file.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.fchk\n";
        return 1;
    }
    std::ifstream from(args[1]);
    if (!from) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }

    const std::string pattern = "Cartesian Force Constants";

    bool found = false;
    double fc;
    std::string line;
    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            found = true;
            int count = 0;
            std::cout << " HESSIAN\n";
            while (from >> fc) {
                std::cout << std::setw(16) << std::setprecision(8)
                          << std::setiosflags(std::ios_base::scientific |
                                              std::ios_base::uppercase)
                          << fc;
                count++;
                if (count == 5) {
                    std::cout << '\n';
                    count = 0;
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

