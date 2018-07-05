////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013-2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/troe.h>
#include <chem/whitten_rabino.h>
#include <srs/datum.h>
#include <srs/utils.h>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <map>

Troe::Troe(std::istream& from, Molecule& mol_) : mol(mol_)
{
    int pot_type_tmp;

    std::map<std::string, srs::Input> input_data;
    input_data["pot_type"]    = srs::Input(pot_type_tmp, 1);
    input_data["e_barrier"]   = srs::Input(e_barrier);
    input_data["imom_ratio"]  = srs::Input(imom_ratio, 1.0);
    input_data["n_free_rot"]  = srs::Input(n_free_rot, 0);
    input_data["n_morse_osc"] = srs::Input(n_morse_osc, 0);

    if (srs::find_section(from, "Troe")) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }
    else {
        throw Troe_error("could not find Troe section");
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Troe_error(it->first + " not initialized");
        }
    }

    // Check if data are sensible:

    if (pot_type_tmp == 1) {
        pot_type = type1;
    }
    else if (pot_type_tmp == 2) {
        pot_type = type2;
    }
    else {
        throw Troe_error("bad potential type");
    }
    if (e_barrier <= 0.0) {
        throw Troe_error("bad energy barrier");
    }
    if (imom_ratio < 0.0) {
        throw Troe_error("bad moment of inertia ratio");
    }
    if (n_morse_osc < 0) {
        throw Troe_error("bad number of Morse oscillators");
    }

    zpe = mol.get_vib().zero_point_energy();
}

double Troe::f_energy(const double temp) const
{
    using namespace boost::math;

    double en = e_barrier + wr::a_corr(mol, e_barrier) * zpe;

    en = datum::R * 1.0e-3 * temp / (en * datum::icm_to_kJ);

    auto s     = mol.get_vib().get_freqs().size();
    double f_e = 0.0;

    for (srs::size_t i = 0; i < s; ++i) {
        f_e += factorial<double>(s - 1) * std::pow(en, i)
               / factorial<double>(s - 1 - i);
    }
    return f_e;
}

double Troe::f_rotation(const double temp) const
{
    using namespace boost::math;
    const double kT = datum::R * 1.0e-3 * temp / datum::icm_to_kJ;

    double e0_azpe = e_barrier + wr::a_corr(mol, e_barrier) * zpe;
    double e0_kT   = e_barrier / kT;
    double f_rot   = 1.0;
    unsigned s     = mol.get_vib().get_freqs().size();

    switch (pot_type) {
    case type1:  // no barrier
        if (mol.structure() == linear) {
            e0_azpe /= s * kT;
            e0_kT = 2.15 * std::pow(e0_kT, 1. / 3.);
            f_rot = e0_azpe * (e0_kT / (e0_kT - 1.0 + e0_azpe));
        }
        else {  // nonlinear
            e0_kT = std::pow(e0_kT, 1. / 3.);
            f_rot = (factorial<double>(s - 1) / tgamma<double>(s + 0.5 + 1.0))
                    * std::pow(e0_azpe / kT, 1.5)
                    * (2.15 * e0_kT
                       / (2.15 * e0_kT - 1.0 + e0_azpe / ((s + 0.5) * kT)));
        }
        break;
    case type2:  // barrier
        if (mol.structure() == linear) {
            e0_azpe /= s * kT;
            f_rot = e0_azpe * imom_ratio / (imom_ratio - 1.0 - e0_azpe);
        }
        else {  // nonlinear
            f_rot = (factorial<double>(s - 1) / tgamma<double>(s + 0.5 + 1.0))
                    * std::pow(e0_azpe / kT, 1.5)
                    * (imom_ratio
                       / (imom_ratio - 1.0 + e0_azpe / ((s + 0.5) * kT)));
        }
        break;
    }
    return f_rot;
}

double Troe::f_free_rotor(const double temp) const
{
    using namespace boost::math;
    double f_free_rot = 1.0;

    if (n_free_rot > 0) {
        double kT = datum::R * 1.0e-3 * temp / datum::icm_to_kJ;

        double e0_azpe = (e_barrier + wr::a_corr(mol, e_barrier) * zpe) / kT;

        auto s = mol.get_vib().get_freqs().size();
        auto r = n_free_rot;

        f_free_rot = (factorial<double>(s - 1) / tgamma<double>(s + 0.5 * r))
                     * std::pow(e0_azpe, 0.5 * r);
    }
    return f_free_rot;
}
