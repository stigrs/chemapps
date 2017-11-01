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

#include <srs/array.h>
#include <srs/math.h>
#include <iostream>
#include <stdexcept>
#include <string>

//-----------------------------------------------------------------------------

// Error reporting:

struct Molvib_error : std::runtime_error {
    Molvib_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

/// Class for handling molecular vibrations.
class Molvib {
public:
    Molvib() {}
    Molvib(std::istream& from, const std::string& key);

    Molvib(const Molvib& vib) { freqs = vib.freqs; }

    /// Get vibrational frequencies.
    const srs::dvector& get_freqs() const { return freqs; }

    /// Calculate zero-point vibrational energy.
    double zero_point_energy() const { return 0.5 * srs::sum(freqs); }

    /// Print vibrational frequencies.
    void print(std::ostream& to = std::cout);

private:
    srs::dvector freqs;
};

#endif  // CHEM_MOLVIB_H
