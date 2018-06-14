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

#include <chem/molvib.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <cmath>


void Molvib::print(std::ostream& to)
{
    srs::Format<char> line;
    line.width(26).fill('-');

    srs::Format<double> fix;
    fix.fixed().width(8).precision(2);

    if (freqs.size() > 0) {
        int it = 0;
        to << "\nVibrational modes (cm^-1):\n" << line('-') << '\n';
        for (srs::size_t i = 0; i < freqs.size(); ++i) {
            to << fix(freqs(i));
            if ((it == 8) && (freqs.size() > 9)) {
                to << '\n';
            }
            it += 1;
        }
        double zpe = zero_point_energy();
        to << "\n\nZero-point vibrational energy: " << zpe / datum::au_to_icm
           << " Hartree\n";
    }
}

srs::packed_dmatrix Molvib::get_mw_hessians() const
{
    srs::packed_dmatrix hess_mw(hess);
    if (!hess_mw.empty()) {
        for (srs::size_t j = 0; j < hess_mw.cols(); ++j) {
            for (srs::size_t i = 0; i < j + 1; ++i) {
                hess_mw(i, j) /= std::sqrt(atoms[i / 3].atomic_mass
                                           * atoms[j / 3].atomic_mass);
            }
        }
    }
    return hess_mw;
}

srs::dvector Molvib::calc_cart_freqs() const
{
    using namespace datum;

    const double factor  // conversion from atomic units to cm-1
        = 0.1 * (N_A * E_h / (4.0 * std::pow(pi * c_0 * a_0 * 1.0e-10, 2.0)));

    srs::packed_dmatrix hess_mw = get_mw_hessians();

    srs::dvector cart_freqs;
    srs::dmatrix v;

    srs::eigs(hess_mw, v, cart_freqs);

    for (auto& vi : cart_freqs) {
        auto x = vi * factor;
        vi     = srs::sign(std::sqrt(std::abs(x)), x);
    }
    return cart_freqs;
}

void Molvib::init(std::istream& from, const std::string& key)
{
    bool found = srs::find_section(from, key);
    if (found) {
        std::string token;
        srs::dvector tmp;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "vibr") {
                from >> freqs;
            }
            else if (token == "hessians") {
                from >> tmp;
                hess = srs::packed_dmatrix(tmp);
            }
        }
    }
    else {
        throw Molvib_error("cannot find " + key + " section");
    }
    // TODO (stigrs@gmail.com) Validate input data
}
