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

#include <fstream>
#include <gsl/gsl>
#include <iostream>


//
// Program for calculation of effective N(E,J) = N1 * N2 / (N1 + N2).
//
int main(int argc, char* argv[])
{
    auto args = gsl::multi_span<char*>(argv, argc);
    if (argc != 3) {
        std::cerr << "usage: " << args[0] << " flux_file_1 flux_file_2\n";
        return 1;
    }

    std::ifstream from1(args[1]);
    if (!from1) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }
    std::ifstream from2(args[2]);
    if (!from2) {
        std::cerr << "cannot open " << args[2] << '\n';
        return 1;
    }

    double e1;
    double e2;
    double j1;
    double j2;
    double n1;
    double n2;
    double neff;

    while (from1 >> e1 >> j1 >> n1) {
        from2 >> e2 >> j2 >> n2;
        if (!from2) {
            std::cerr << "input flux file " << args[2] << " too short?\n";
            return 1;
        }
        if (e1 != e2) {
            std::cerr << "bad E grid: " << e1 << ", " << e2 << '\n';
            return 1;
        }
        if (j1 != j2) {
            std::cerr << "bad J grid: " << j1 << ", " << j2 << '\n';
            return 1;
        }
        if ((n1 <= 0.0) || (n2 <= 0.0)) {
            neff = 0.0;
        }
        else {
            neff = n1 * n2 / (n1 + n2);
        }
        std::cout << e1 << " " << j1 << " " << neff << '\n';
    }
    if (from2 >> n2) {
        std::cerr << "input flux file " << args[1] << " too short?\n";
        return 1;
    }
}
