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

#ifndef CHEM_MOLVIB_H
#define CHEM_MOLVIB_H

#include <chem/element.h>
#include <srs/array.h>
#include <srs/math.h>
#include <srs/packed.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------

// Error reporting:

struct Molvib_error : std::runtime_error {
    Molvib_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

// Class for handling molecular vibrations.
class Molvib {
public:
    Molvib() {}
    Molvib(std::istream& from,
           const std::string& key,
           std::vector<Element>& atoms_)
        : atoms(atoms_)
    {
        init(from, key);
    }

    Molvib(const Molvib& vib)
        : atoms(vib.atoms), freqs(vib.freqs), hess(vib.hess)
    {
    }

    // Get vibrational frequencies.
    const srs::dvector& get_freqs() const { return freqs; }

    // Get Hessians.
    const srs::packed_dmatrix& get_hessians() const { return hess; }

	// Get mass-weighted Hessians.
    srs::packed_dmatrix get_mw_hessians() const;

	// Calculate Cartesian frequencies from Hessians.
    srs::dvector calc_cart_freqs() const;

    // Calculate zero-point vibrational energy.
    double zero_point_energy() const { return 0.5 * srs::sum(freqs); }

    // Print vibrational frequencies.
    void print(std::ostream& to = std::cout);

private:
    // Initialize.
    void init(std::istream& from, const std::string& key);

    std::vector<Element> atoms;
    srs::dvector freqs;
    srs::packed_dmatrix hess;
};

#endif  // CHEM_MOLVIB_H
