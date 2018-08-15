////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/energy_levels.h>
#include <srs/math.h>
#include <cmath>
#include <gsl/gsl>

srs::dvector energy_levels::harmonic_oscillator(double freq, double emax)
{
    Expects(freq > 0.0);
    Expects(emax > 0.0);

    int kmax = 1 + srs::round<int>(emax / freq);
    srs::dvector result(kmax);
    for (int i = 0; i < kmax; ++i) {
        result(i) = freq * (i + 1);
    }
    return result;
}

srs::dvector energy_levels::free_rotor(double rotc, double emax)
{
    Expects(rotc > 0.0);
    Expects(emax > 0.0);

    srs::dvector result;
    double ej = 0.0;
    int j     = 1;
    while (ej < emax) {
        ej = rotc * j * j;
        result.push_back(ej);
        ++j;
    }
    return result;
}

srs::dvector energy_levels::hindered_rotor(double sigma,
                                           double rotc,
                                           double barrier,
                                           double emax)
{
    Expects(sigma >= 1.0);
    Expects(rotc > 0.0);
    Expects(barrier >= 0.0);
    Expects(emax > 0.0);

    srs::dvector result;

    if (barrier > 1.0) {  // hindered rotor
        double sig = sigma;
        double b   = rotc;
        double v0  = barrier;
        double frq = sig * std::sqrt(b * v0);
        double r   = v0 / frq;
        double es  = 0.0;
        double zpe = 0.0;
        int ns     = 0;
        while (es < emax) {
            // Calculate harmonic oscillator energy level:
            double dnv = ns / sig;
            int nv     = srs::round<int>(dnv);
            double tv  = -frq * (1.0 + 2.0 * nv + 2.0 * nv * nv) / (16.0 * r);
            double ev  = frq * (nv + 0.5) + tv;

            // Calculate free rotor energy level:
            int j     = (ns + 1) / 2;
            double tr = 0.0;
            if ((j > (r * sig / 2.0)) && (r > 0.0)) {
                tr = std::pow(r, 4.0) * sig * sig * b
                     / (8.0 * (std::pow(2.0 * j / sig, 2.0) - 1.0));
            }
            double ej = b * j * j + 0.5 * v0 + tr;

            // Calculate hindered rotor energy level:
            double s = 0.5 * (1.0 + std::tanh(5.0 * (ev - v0) / v0));
            if (ej > 1.5 * v0) {
                s = 1.0;
            }
            es = ev * (1.0 - s) + ej * s;

            if (ns == 0) {
                zpe = es;
            }
            if (ns > 0) {
                result.push_back(es - zpe);
            }
            ++ns;
        }
    }
    else {  // free rotor
        result = free_rotor(rotc, emax);
    }
    return result;
}
