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

#include <chem/tunnel.h>
#include <numlib/math.h>
#include <array>
#include <string>
#include <stdexcept>

Chem::Tunnel::Tunnel(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    std::string method_str;
    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "method", method_str, std::string("none"));
        get_token_value(from, pos, "freq_im", freq_im, 0.0);
        get_token_value(from, pos, "en_barrier", en_barrier, 0.0);
        get_token_value(from, pos, "en_rxn", en_rxn, 0.0);
    }
    Assert::dynamic(freq_im < 0.0, "bad freq_im");

    // Set tunneling correction method:
    if (method_str == "None") {
        method = None;
    }
    else if (method_str == "Wigner") {
        method = Wigner;
    }
    else if (method_str == "Eckart") {
        method = Eckart;
    }
    else {
        throw std::runtime_error("unknown tunneling correction: " + method_str);
    }
}

double Chem::Tunnel::eckart(double temp) const
{
    // The implementation is based on the following papers:
    //
    //  1. Eckart, E. Phys. Rev., 1962, vol. 35, p. 1303.
    //  2. Brown, R. L. J. Research NIST, 1981, vol. 86, p. 357.
    //  3. Johnston, H. S.; Heicklen, J. J. Phys. Chem., 1962, vol. 66, p. 532.
    //
    // Algorithm:
    // ----------
    // Brown (1981) introduced a new variable
    //
    //     epsilon = (energy - en_barrier) / (k_boltzmann * temp)
    //
    // in order to evaluate the integral yielding the tunneling correction.
    // When epsilon gets large, the transmission probability, kappa, approaches
    // unity. This method uses a Gaussian formula for the part of the integral
    // where kappa < 1. The remainder where kappa is ca. 1 is evaluated
    // analytically. The energy at which this happens is called epsilon_b and
    // kappa is called kappa_b. However, epsilon_b is kept below a certain
    // constant value epsilon_max. The parameters kappa_b and epsilon_max have
    // been adjusted to give best possible agreement with the values presented
    // in Table I in the paper of Johnston and Heicklen (1962). The difference
    // is less than 1 percent and therefore should give better results than
    // the program of Brown (1962). This has not been fully tested though.

    using namespace Numlib::Constants;

    Assert::dynamic<Assert::level(2)>(temp > 0.0, "bad temperature");

    double kt = k * temp;
    double ifreq = std::abs(freq_im) * c_0 * 100.0;

    double pot1 = en_barrier * kilo / N_A;
    double pot2 = (en_barrier - en_rxn) * kilo / N_A;

    double alpha1 = 2.0 * pi * pot1 / (h * ifreq);
    double alpha2 = 2.0 * pi * pot2 / (h * ifreq);

    double u = h * ifreq / kt;
    double d = 4.0 * alpha1 * alpha2 - pi * pi;
    double c = 0.125 * pi * u *
               std::pow(1.0 / std::sqrt(alpha1) + 1.0 / std::sqrt(alpha2), 2.0);

    double df = std::cosh(std::sqrt(d));
    if (d <= 0.0) {
        df = std::cos(std::sqrt(std::abs(d)));
    }

    constexpr double kappa_b = 1.0e-10;   // this is actually 1 - kappa_b
    constexpr double epsilon_max = 110.0; // not the same as Brown (1981)

    double epsilon_b =
        c * std::pow(std::log(2.0 * (1.0 + df) / kappa_b) / (2.0 * pi), 2.0);
    epsilon_b = std::min(epsilon_b, epsilon_max);

    double v1 = pot1 / kt;
    double v2 = pot2 / kt;

    double epsilon_0 = 0.0;
    if (v2 >= v1) {
        epsilon_0 = -v1;
    }
    if (v1 > v2) {
        epsilon_0 = -v2;
    }

    // Integrate using a 16-point Gauss-Legendre quadrature:

    const int n = 16;
    std::array<double, n> x;
    std::array<double, n> w;
    Numlib::gauss_legendre<n>(x, w, epsilon_0, epsilon_b);

    double kappa = 0.0;
    for (int i = 0; i < n; ++i) {
        double a1 = pi * std::sqrt((x[i] + v1) / c);
        double a2 = pi * std::sqrt((x[i] + v2) / c);
        double fp = std::cosh(a1 + a2);
        double fm = std::cosh(a1 - a2);
        kappa += w[i] * std::exp(-x[i]) * (fp - fm) / (fp + df);
    }

    // Add the analytic part:
    kappa += std::exp(-epsilon_b);
    return kappa;
}

