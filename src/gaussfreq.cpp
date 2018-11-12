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

#include <chem/gauss_data.h>
#include <stdutils/stdutils.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

//  Extracts vibrational frequencies from Gaussian output file.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " gaussian.log\n";
        return 1;
    }

    std::ifstream from(args[1]);
    if (!from) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }

    std::vector<double> freqs;

    Chem::Gauss_data gauss(from, Chem::out);
    gauss.get_freqs(freqs);

    for (auto vi : freqs) {
        std::cout << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4) << vi << '\t';
    }
    std::cout << '\n';
}
