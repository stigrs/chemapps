//////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
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
///////////////////////////////////////////////////////////////////////////////

#ifndef CHEM_ARMA_IO_H
#define CHEM_ARMA_IO_H

#include <armadillo>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace chem {

// Error reporting:

struct Arma_error : std::runtime_error {
    Arma_error(std::string s) : std::runtime_error(s) {}
};

// I/O operators for vectors:

template <typename T>
void print_vector(std::ostream& to, const arma::Col<T>& a)
{
    to << a.size() << " [ ";
    for (arma::uword i = 0; i < a.size(); ++i) {
        to << a(i) << " ";
        if (!((i + 1) % 7) && (i != (a.size() - 1))) {
            to << "\n  ";
        }
    }
    to << ']';
}

template <typename T>
void read_vector(std::istream& from, arma::Col<T>& a)
{
    arma::uword n;
    from >> n;
    if (n < 1) {
        throw Arma_error("read_vector: bad size");
    }
    a.set_size(n);

    char ch;
    from >> ch;
    if (ch != '[') {
        throw Arma_error("read_vector: '[' missing");
    }
    for (arma::uword i = 0; i < n; ++i) {
        from >> a(i);
    }
    from >> ch;
    if (ch != ']') {
        throw Arma_error("read_vector: ']' missing");
    }
}
}  // namespace chem

#endif  // CHEM_ARMA_IO_H
