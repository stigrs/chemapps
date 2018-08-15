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
#include <chem/statecount.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>


srs::dvector statecount::count(const Molecule& mol,
                               int ngrains,
                               double egrain,
                               bool sum)
{
    srs::dvector rot;
    if (mol.has_torsions()) {
        double sigma   = mol.get_tor().symmetry_number();
        double rotc    = srs::max(mol.get_tor().constant());
        double barrier = srs::max(mol.get_tor().get_pot_coeff());
        if (barrier < 0.01) {  // free rotor
            rot = free_rotor(sigma, rotc, ngrains, egrain, sum);
        }
        else {  // hindered rotor
            rot = hindered_rotor(sigma, rotc, barrier, ngrains, egrain, sum);
        }
    }
    return bswine(mol.get_vib().get_freqs(), ngrains, egrain, sum, rot);
}

srs::dvector statecount::bswine(const srs::dvector& vibr,
                                int ngrains,
                                double egrain,
                                bool sum,
                                const srs::dvector& rot)
{
    srs::dvector result = srs::zeros<srs::dvector>(ngrains);
    if (!rot.empty()) {  // initialize with rotational states
        result = rot;
    }
    else {
        if (sum) {  // count sum of states
            result = srs::ones<srs::dvector>(ngrains);
        }
        else {  // count density of states
            result(0) = 1.0;
        }
    }
    for (auto w : vibr) {
        int wj = srs::round<int>(w / egrain);
        for (int i = wj; i < ngrains; ++i) {
            result(i) += result(i - wj);
        }
    }
    if (!sum) {
        result *= 1.0 / egrain;
    }
    return result;
}

srs::dvector statecount::steinrab(const srs::dvector& vibr,
                                  double sigma,
                                  double rotc,
                                  double barrier,
                                  int ngrains,
                                  double egrain,
                                  bool sum)
{
    srs::dvector at = srs::zeros<srs::dvector>(ngrains);
    srs::dvector tt = srs::zeros<srs::dvector>(ngrains);

    at(0) = 1.0;
    tt(0) = 1.0;

    double emax = ngrains * egrain;

    if (rotc != 0.0) {
        srs::dvector rj;
        double dd = 0.0;      // degeneracy
        if (barrier > 1.0) {  // hindered rotor
            dd = 1.0;
            rj = energy_levels::hindered_rotor(sigma, rotc, barrier, emax);
        }
        else {  // free rotor
            dd = 2.0;
            rj = energy_levels::free_rotor(rotc, emax);
        }
        for (srs::size_t k = 0; k < rj.size(); ++k) {
            int rjk = srs::round<int>(rj(k) / egrain);
            for (int i = rjk; i < ngrains; ++i) {
                at(i) += dd * tt(i - rjk);
            }
        }
        tt = (1.0 / sigma) * at;
        at = tt;
    }
    for (auto w : vibr) {
        auto rj = energy_levels::harmonic_oscillator(w, emax);
        for (srs::size_t k = 0; k < rj.size(); ++k) {
            int rjk = srs::round<int>(rj(k) / egrain);
            for (int i = rjk; i < ngrains; ++i) {
                at(i) += tt(i - rjk);
            }
        }
        tt = at;
    }
    if (sum) {
        for (int i = 1; i < ngrains; ++i) {
            tt(i) += tt(i - 1);
        }
    }
    else {
        tt *= 1.0 / egrain;
    }
    return tt;
}

srs::dvector statecount::free_rotor(
    double sigma, double rotc, int ngrains, double egrain, bool sum)
{
    srs::dvector result = srs::zeros<srs::dvector>(ngrains);

    double f = -0.5;
    if (sum) {
        f = 0.5;
    }
    double qr = std::sqrt(datum::pi) / (sigma * std::sqrt(rotc));
    double gf = boost::math::tgamma(f + 1.0);
    for (int i = 0; i < ngrains; ++i) {
        result(i) = qr * std::pow(i * egrain, f) / gf;
    }
    return result;
}

srs::dvector statecount::hindered_rotor(double sigma,
                                        double rotc,
                                        double barrier,
                                        int ngrains,
                                        double egrain,
                                        bool sum)
{
    using namespace boost::math;

    srs::dvector result = srs::zeros<srs::dvector>(ngrains);

    int iv0 = srs::round<int>(barrier / egrain - 0.5);

    double v0  = barrier;
    double q1f = std::sqrt(datum::pi) / (sigma * std::sqrt(rotc));
    double ei  = 0.0;

    if (sum) {  // eq. 4.52 in Forst (2003)
        for (int i = 0; i < iv0; ++i) {
            ei = i * egrain;
            result(i)
                = (4.0 * q1f * std::sqrt(v0) / std::pow(datum::pi, 1.5))
                  * (ellint_2(ei / v0) - (1.0 - ei / v0) * ellint_1(ei / v0));
        }
        // Ugly hack for E/V0 = 1 in order to avoid NaN:
        ei          = iv0 * egrain;
        result(iv0) = (4.0 * q1f * std::sqrt(v0) / std::pow(datum::pi, 1.5))
                      * ellint_2(ei / v0);
        for (int i = iv0 + 1; i < ngrains; ++i) {
            ei        = i * egrain;
            result(i) = (4.0 * q1f * std::sqrt(ei) / std::pow(datum::pi, 1.5))
                        * ellint_2(v0 / ei);
        }
    }
    else {
        for (int i = 0; i < iv0; ++i) {
            ei        = i * egrain;
            result(i) = (2.0 * q1f / (std::pow(datum::pi, 1.5) * std::sqrt(v0)))
                        * ellint_1(ei / v0);
        }
        for (int i = iv0 + 1; i < ngrains; ++i) {
            ei        = i * egrain;
            result(i) = (2.0 * q1f / (std::pow(datum::pi, 1.5) * std::sqrt(ei)))
                        * ellint_1(v0 / ei);
        }
        // Ugly hack for E/V0 = 1 in order to avoid Inf:
        result(iv0) = 0.5 * (result(iv0 + 1) + result(iv0 - 1));
    }
    return result;
}
