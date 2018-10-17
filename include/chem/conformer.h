////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

#ifndef CHEM_CONFORMER_H
#define CHEM_CONFORMER_H

#include <chem/element.h>
#include <numlib/matrix.h>
#include <vector>

namespace Chem {

// Simple structure for storing conformers.
//
struct Conformer {
    Conformer() = default;

    explicit Conformer(double e, const Numlib::Mat<double>& x)
        : energy{e}, atoms{0} xyz(x)
    {
        iter = 0;
    }

    Conformer(double e,
              const std::vector<Element>& at,
              const Numlib::Mat<double>& x)
        : energy{e}, atoms(at), xyz(x)
    {
        iter = 0;
    }

    // Copy semantics:
    Conformer(const Conformer&) = default;
    Conformer& operator=(const Conformer&) = default;

    // Move semantics:
    Conformer(Conformer&&) = default;
    Conformer& operator=(Conformer&&) = default;

    ~Conformer() = default;

    // Compare conformers by energy.
    bool operator<(const Conformer& c) const { return energy < c.energy; }

    // Compare conformers by energy.
    bool operator>(const Conformer& c) const { return energy > c.energy; }

    double energy;               // conformer energy
    std::vector<Element> atoms;  // atoms
    Numlib::Mat<double> xyz;     // Cartesian coordinates
    int iter;                    // iterator to be used for counting
};

} // namespace Chem

#endif // CHEM_CONFORMER_H

