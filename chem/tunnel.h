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

#ifndef CHEM_TUNNEL_H
#define CHEM_TUNNEL_H

#include <chem/datum.h>
#include <cmath>
#include <gsl/gsl>
#include <iostream>
#include <stdexcept>
#include <string>

// Error reporting:

struct Tunnel_error : std::runtime_error {
    Tunnel_error(std::string s) : std::runtime_error(s) {}
};

//
// Class for computing quantum tunneling corrections.
//
class Tunnel {
public:
    Tunnel(std::istream& from, const std::string& key = "Tunnel");

    ~Tunnel() {}

    // Calculate Wigner tunneling correction.
    double wigner(double temp = 298.15) const;

    // Calculate Eckart tunneling correction for an unsymmetrical barrier.
    double eckart(double temp = 298.15) const;

    // Calculate tunneling correction factor.
    double factor(double temp = 298.15) const;

private:
    enum Method_t { Wigner, Eckart };

    Method_t method = Wigner;  // tunneling correction method
    double freq_im;            // imaginary frequency
    double en_barrier;         // potential barrier height
    double en_rxn;             // energy of reaction
};

inline double Tunnel::wigner(double temp) const
{
    // Wigner, E. Z. Physik. Chem. (Leipzig), 1932, vol. B19, p. 203.

    Expects(temp > 0.0);
    double factor
        = datum::h * std::abs(freq_im) * 100.0 * datum::c_0 / (datum::k * temp);
    return 1.0 + factor * factor / 24.0;
}

inline double Tunnel::factor(double temp) const
{
    switch (method) {
    case Wigner:
        return wigner(temp);
    case Eckart:
        return eckart(temp);
    default:
        return 1.0;
    }
}

#endif  // CHEM_TUNNEL_H
