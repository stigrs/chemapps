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

#ifndef CHEM_CONFORMER_H
#define CHEM_CONFORMER_H

#include <chem/element.h>
#include <armadillo>
#include <vector>

//
// Simple structure for storing conformers.
//
struct Conformer {
    Conformer() {}
    explicit Conformer(double e, const arma::mat& x)
        : energy(e), atoms(0), xyz(x)
    {
        iter = 0;
    }

    Conformer(double e, const std::vector<Element>& at, const arma::mat& x)
        : energy(e), atoms(at), xyz(x)
    {
        iter = 0;
    }

    Conformer(const Conformer& c) : atoms(c.atoms), xyz(c.xyz)
    {
        energy = c.energy;
        iter   = c.iter;
    }

    // Compare conformers by energy.
    bool operator<(const Conformer& c) const { return energy < c.energy; }

    // Compare conformers by energy.
    bool operator>(const Conformer& c) const { return energy > c.energy; }

    double energy;               // conformer energy
    std::vector<Element> atoms;  // atoms
    arma::mat xyz;               // Cartesian coordinates
    int iter;                    // iterator to be used for counting
};

#endif  // CHEM_CONFORMER_H
