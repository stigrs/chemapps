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
#include <chem/utils.h>
#include <string>

void chem::thermochemistry(const Molecule& mol,
                           const arma::vec& temp,
                           const arma::vec& pressure,
                           bool incl_sigma,
                           std::ostream& to)
{
    chem::Format<char> line;
    line.width(16).fill('=');

    chem::Format<double> fix;
    fix.fixed().precision(6);

    to << "\nThermochemistry:\n" << line('=') << '\n';
    to << "Electronic energy: " << fix(mol.get_elec_energy()) << " Hartree\n";

    if (mol.has_torsions()) {
        mol.get_tor().analysis(to);
    }
    mol.get_vib().print(to);
    to << '\n';

    fix.fixed().precision(3);

    typedef arma::vec::const_iterator Citer;  // arma:: does not support auto
    for (Citer p = pressure.begin(); p != pressure.end(); ++p) {
        for (Citer t = temp.begin(); t != temp.end(); ++t) {
            to << "Temperature: " << fix(*t) << " K. "
               << "Pressure: " << fix(*p) << " bar\n";
        }
    }
    incl_sigma;
}

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

    if (rot_symm.find("atom") != std::string::npos) {
        return 1.0;
    }
    else if (rot_symm.find("linear") != std::string::npos) {
        arma::vec3 rotc = datum::GHz_to_K * mol.get_rot().constants();
        double rsig     = 1.0;
        if (incl_sigma) {
            rsig /= mol.get_rot().get_sigma();
        }
        return rsig * temp / rotc(0);
    }
    else {  // nonlinear molecule
        arma::vec3 rotc = datum::GHz_to_K * mol.get_rot().constants();
        double b        = arma::prod(rotc);
        double rsig     = std::sqrt(datum::pi);
        if (incl_sigma) {
            rsig /= mol.get_rot().get_sigma();
        }
        return rsig * std::pow(temp, 1.5) / std::sqrt(b);
    }
}

double chem::entropy_rot(const Molecule& mol, double temp, bool incl_sigma)
{
    std::string rot_symm = mol.get_rot().symmetry();

    if (rot_symm.find("atom") != std::string::npos) {
        return 1.0;
    }
    else {
        double factor = 1.5;
        if (rot_symm.find("linear") != std::string::npos) {
            factor = 1.0;
        }
        double qr = chem::qrot(mol, temp, incl_sigma);
        Ensures(qr > 0.0);
        return datum::R * (std::log(qr) + factor);
    }
}

double chem::qvib(const Molecule& mol, double temp, const std::string& zeroref)
{
    std::string rot_symm = mol.get_rot().symmetry();

    if (rot_symm.find("atom") != std::string::npos) {
        return 1.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        double qv   = 1.0;
        if (zeroref == "V=0") {  // zero at first vibrational level
            for (arma::uword i = 0; i < w.size(); ++i) {
                qv /= (1.0 - std::exp(-w(i) / temp));
            }
        }
        else if (zeroref == "BOT") {  // zero at the bottom of the well
            for (arma::uword i = 0; i < w.size(); ++i) {
                qv *= std::exp(-w(i) / (2.0 * temp))
                      / (1.0 - std::exp(-w(i) / temp));
            }
        }
        return qv;
    }
}

double chem::entropy_vib(const Molecule& mol, double temp)
{
    std::string rot_symm = mol.get_rot().symmetry();

    if (rot_symm.find("atom") != std::string::npos) {
        return 0.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        double sv   = 0.0;
        for (arma::uword i = 0; i < w.size(); ++i) {
            double wt = w(i) / temp;
            sv += wt / (std::exp(wt) - 1.0) - std::log(1.0 - std::exp(-wt));
        }
        sv *= datum::R;
        return sv;
    }
}

double chem::thermal_energy_vib(const Molecule& mol, double temp)
{
    std::string rot_symm = mol.get_rot().symmetry();

    if (rot_symm.find("atom") != std::string::npos) {
        return 0.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        double ev   = 0.0;
        for (arma::uword i = 0; i < w.size(); ++i) {
            ev += w(i) * (0.5 + 1.0 / (std::exp(w(i) / temp) - 1.0));
        }
        ev *= datum::R;
        return ev;
    }
}

double chem::const_vol_heat_vib(const Molecule& mol, double temp)
{
    std::string rot_symm = mol.get_rot().symmetry();

    if (rot_symm.find("atom") != std::string::npos) {
        return 0.0;
    }
    else {
        Expects(temp > 0.0);
        arma::vec w = datum::icm_to_K * mol.get_vib().get_freqs();
        double cv_v = 0.0;
        for (arma::uword i = 0; i < w.size(); ++i) {
            double wt = w(i) / temp;
            cv_v += wt * wt * std::exp(wt) / std::pow(std::exp(wt) - 1.0, 2.0);
        }
        cv_v *= datum::R;
        return cv_v;
    }
}

double chem::qctcw(const Molecule& mol, double temp)
{
    double qtor = 1.0;

    std::string rot_symm = mol.get_rot().symmetry();
    if (rot_symm.find("atom") != std::string::npos) {
        qtor = 1.0;
    }
    else {
        if (mol.has_torsions()) {
            Expects(temp > 0.0);
            // Calculate free rotor partition function:
            double imom = mol.get_tor().eff_moment_of_inertia();
            imom *= datum::au_to_kgm2;
            double sigma = mol.get_tor().symmetry_number();
            double qfr   = std::sqrt(2.0 * datum::pi * imom * datum::k * temp)
                         / (datum::h_bar * sigma);

            // Calculate partition function for harmonic oscillator and
            // intermediate case:
            double qho     = 0.0;
            double qin     = 0.0;
            arma::vec pot  = mol.get_tor().get_pot_coeff();
            arma::vec freq = mol.get_tor().get_freqs();
            Expects(pot.size() == freq.size());
            for (arma::uword i = 0; i < pot.size(); ++i) {
                double ui = pot(i) * datum::icm_to_K;
                double wi = freq(i) * datum::icm_to_K;
                qho += std::exp(-(ui + 0.5 * wi) / temp)
                       / (1.0 - std::exp(-wi / temp));
                qin += std::exp(-ui / temp) / (wi / temp);
            }
            qtor = qho * std::tanh(qfr / qin);  // eq. 11 in C&T (2000).
        }
    }
    return qtor;
}

double chem::const_vol_heat_tor(const Molecule& mol, double temp)
{
    // The constant volume heat capacity is calculated as the numerical
    // derivative of the thermal torsional energy with respect to
    // temperature (dEtor/dT) at constant N and V.

    std::string rot_symm = mol.get_rot().symmetry();
    if (rot_symm.find("atom") != std::string::npos) {
        return 0.0;
    }
    else {
        if (rot_symm.find("linear") != std::string::npos) {
            return 0.0;  // a linear molecule cannot have torsional modes
        }
        else {
            Expects(temp > 0.0);
            double h = temp * std::pow(std::numeric_limits<double>::epsilon(),
                                       1.0 / 3.0);
            double cva = chem::thermal_energy_tor(mol, temp + h);
            double cvb = chem::thermal_energy_tor(mol, temp - h);
            return (cva - cvb) / (2.0 * h);
        }
    }
}