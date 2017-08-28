///////////////////////////////////////////////////////////////////////////////
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

#include <chem/thermochem.h>
#include <armadillo>
#include <string>

double chem::qelec(const Molecule& mol, double temp)
{
    arma::vec elec = mol.get_elec_state();
    double qe      = 0.0;
    for (arma::uword i = 0; i < elec.size(); i += 2) {
        qe += elec(i) * std::exp(-elec(i + 1) * datum::icm_to_K / temp);
    }
    return qe;
}

double chem::qrot(const Molecule& mol, double temp, bool incl_sigma)
{
    Expects(temp >= 0.0);

    std::string rot_symm = mol.get_rot().symmetry();

    double qr = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        qr = 1.0;
    }
    else if (rot_symm.find("linear") != std::string::npos) {
        arma::vec3 rotc = datum::GHz_to_K * mol.get_rot().constants();
        double rsig     = 1.0;
        if (incl_sigma) {
            rsig /= mol.get_rot().get_sigma();
        }
        qr = rsig * temp / rotc(0);
    }
    else {  // nonlinear molecule
        arma::vec3 rotc = datum::GHz_to_K * mol.get_rot().constants();
        double b        = arma::prod(rotc);
        double rsig     = std::sqrt(datum::pi);
        if (incl_sigma) {
            rsig /= mol.get_rot().get_sigma();
        }
        qr = rsig * std::pow(temp, 1.5) / std::sqrt(b);
    }
    return qr;
}

double chem::entropy_rot(const Molecule& mol, double temp, bool incl_sigma)
{
    std::string rot_symm = mol.get_rot().symmetry();

    double sr = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        sr = 1.0;
    }
    else {
        double factor = 1.5;
        if (rot_symm.find("linear") != std::string::npos) {
            factor = 1.0;
        }
        double qr = chem::qrot(mol, temp, incl_sigma);
        Ensures(qr > 0.0);
        sr = datum::R * (std::log(qr) + factor);
    }
    return sr;
}

double chem::qvib(const Molecule& mol, double temp, const std::string& zeroref)
{
    std::string rot_symm = mol.get_rot().symmetry();

    double qv = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        qv = 1.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        qv          = 1.0;
        if (zeroref == "V=0") {  // zero at first vibrational level
            for (arma::uword i = 0; i < w.size(); ++i) {
                qv /= (1.0 - std::exp(-w(i) / temp));
            }
        }
        else if (zeroref == "BOT") {  // zero at the bottom of the well
            for (arma::uword i = 0; i < w.size(); ++i) {
                qv *= std::exp(-w(i) / (2.0 * temp)) /
                      (1.0 - std::exp(-w(i) / temp));
            }
        }
    }
    return qv;
}

double chem::entropy_vib(const Molecule& mol, double temp)
{
    std::string rot_symm = mol.get_rot().symmetry();

    double sv = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        sv = 0.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        for (arma::uword i = 0; i < w.size(); ++i) {
            double wt = w(i) / temp;
            sv += wt / (std::exp(wt) - 1.0) - std::log(1.0 - std::exp(-wt));
        }
        sv *= datum::R;
    }
    return sv;
}

double chem::thermal_energy_vib(const Molecule& mol, double temp)
{
    std::string rot_symm = mol.get_rot().symmetry();

    double ev = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        ev = 0.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        for (arma::uword i = 0; i < w.size(); ++i) {
            ev += w(i) * (0.5 + 1.0 / (std::exp(w(i) / temp) - 1.0));
        }
        ev *= datum::R;
    }
    return ev;
}

double chem::const_vol_heat_vib(const Molecule& mol, double temp)
{
    std::string rot_symm = mol.get_rot().symmetry();

    double cv_v = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        cv_v = 0.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        for (arma::uword i = 0; i < w.size(); ++i) {
            double wt = w(i) / temp;
            cv_v += wt * wt * std::exp(wt) / std::pow(std::exp(wt) - 1.0, 2.0);
        }
        cv_v *= datum::R;
    }
    return cv_v;
}

double chem::qctcw(const Molecule& mol, double temp)
{
    double qtor = 1.0;

    std::string rot_symm = mol.get_rot().symmetry();
    if (rot_symm.find("atom") != std::string::npos) {
        qtor = 1.0;
    }
    return qtor * temp;
}