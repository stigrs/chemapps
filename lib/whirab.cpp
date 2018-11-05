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

#include <chem/whirab.h>
#include <chem/traits.h>
#include <numlib/constants.h>
#include <numlib/traits.h>
#include <cmath>

double Chem::Whirab::a_corr(const Chem::Molecule& mol, double e_barrier)
{
    double en = e_barrier / mol.zero_point_energy();
    double w;
    if (en >= 1.0) {
        w = -1.0506 * std::pow(en, 0.25);
        w = std::pow(10.0, w);
    }
    else {
        w = 1.0 / (5.0 * en + 2.73 * std::sqrt(en) + 3.51);
    }
    double sum2_v = 0.0;
    double sum_v2 = 0.0;
    for (auto vi : mol.frequencies()) {
        if (vi < 0.0) {  // ignore imaginary frequencies
            continue;
        }
        else {
            sum2_v += vi;
            sum_v2 += vi * vi;
        }
    }
    sum2_v *= sum2_v;

    double factor = 3.0;  // structure factor for nonlinear
    if (mol.structure() == linear) {
        factor = 2.0;
    }
    auto s = narrow_cast<double>(mol.frequencies().size());
    auto r = narrow_cast<double>(mol.tor_pot_coeff().size());

    double beta
        = (s - 1.0) * ((s + 0.5 * r + 0.5 * factor) / s) * sum_v2 / sum2_v;

    return 1.0 - beta * w;
}

double Chem::Whirab::vibr_density_states(const Chem::Molecule& mol,
                                         double e_barrier)
{
    double hv = 1.0;
    for (auto wi : mol.frequencies()) {
        if (wi < 0.0) {  // ignore imaginary frequencies
            continue;
        }
        else {
            hv *= wi;
        }
    }
    double s = narrow_cast<double>(mol.frequencies().size());
    double sm1 = s - 1.0;
    double rho = e_barrier + a_corr(mol, e_barrier) * mol.zero_point_energy();

    rho = std::pow(rho, sm1) / (std::tgamma(s) * hv);
    return rho / Numlib::Constants::icm_to_kJ;
}

