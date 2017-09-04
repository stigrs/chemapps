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

#include <chem/datum.h>
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

    double e0 = mol.get_elec_energy();
    to << "\nThermochemistry:\n" << line('=') << '\n';
    to << "Electronic energy: " << fix(e0) << " Hartree\n";

    mol.get_rot().analysis(to);

    if (mol.has_torsions()) {
        mol.get_tor().analysis(to);
    }
    mol.get_vib().print(to);
    to << '\n';

    double zpe          = mol.get_vib().zero_point_energy() / datum::au_to_icm;
    const double factor = 1.0 / (datum::E_h * datum::N_A);

    typedef arma::vec::const_iterator Citer;  // arma:: does not support auto
    for (Citer p = pressure.begin(); p != pressure.end(); ++p) {
        for (Citer t = temp.begin(); t != temp.end(); ++t) {
            to << "Temperature: " << *t << " K. "
               << "Pressure: " << *p << " bar\n";

            fix.fixed().width(12).precision(6);
            to << "Zero-point correction:\t\t\t\t" << fix(zpe) << " Hartree\n";

            double ecorr = chem::thermal_energy(mol, *t) * factor;
            to << "Thermal correction to energy:\t\t\t" << fix(ecorr) << '\n';

            double hcorr = chem::enthalpy(mol, *t) * factor;
            to << "Thermal correction to enthalpy:\t\t\t" << fix(hcorr) << '\n';

            double gcorr = chem::gibbs_energy(mol, *t, *p, incl_sigma) * factor;
            to << "Thermal correction to Gibbs energy:\t\t" << fix(gcorr)
               << '\n';

            double etot = e0 + zpe;
            to << "Sum of electronic and zero-point energies:\t" << fix(etot)
               << '\n';

            etot = e0 + ecorr;
            to << "Sum of electronic and thermal energies:\t\t" << fix(etot)
               << '\n';

            double htot = e0 + hcorr;
            to << "Sum of electronic and thermal enthalpies:\t" << fix(htot)
               << '\n';

            double gtot = e0 + gcorr;
            to << "Sum of electronic and Gibbs free energies:\t" << fix(gtot)
               << "\n\n";

            line.width(64).fill('-');
            to << "\t\t\t"
               << "E(thermal)\t"
               << "CV\t\t"
               << "S\t" << '\n'
               << "\t\t\t"
               << "kJ/mol\t\t"
               << "J/mol-K\t\t"
               << "J/mol-K\n"
               << line('-') << '\n';

            fix.width(8).precision(3);
            etot        = ecorr * datum::au_to_icm * datum::icm_to_kJ;
            double cv   = chem::const_vol_heat_capacity(mol, *t);
            double stot = chem::entropy(mol, *t, *p, incl_sigma);
            to << "Total:\t\t\t" << fix(etot) << '\t' << fix(cv) << '\t'
               << fix(stot) << '\n';

            double eelec = chem::thermal_energy_elec();
            double celec = chem::const_vol_heat_elec();
            double selec = chem::entropy_elec(mol, *t);
            to << "Electronic:\t\t" << fix(eelec) << '\t' << fix(celec) << '\t'
               << fix(selec) << '\n';

            double etrans = chem::thermal_energy_trans(*t) / datum::kilo;
            double ctrans = chem::const_vol_heat_trans();
            double strans = chem::entropy_trans(mol, *t, *p);
            to << "Translational:\t\t" << fix(etrans) << '\t' << fix(ctrans)
               << '\t' << fix(strans) << '\n';

            double erot = chem::thermal_energy_rot(mol, *t) / datum::kilo;
            double crot = chem::const_vol_heat_rot(mol);
            double srot = chem::entropy_rot(mol, *t, incl_sigma);
            to << "Rotational:\t\t" << fix(erot) << '\t' << fix(crot) << '\t'
               << fix(srot) << '\n';

            double evib = chem::thermal_energy_vib(mol, *t) / datum::kilo;
            double cvib = chem::const_vol_heat_vib(mol, *t);
            double svib = chem::entropy_vib(mol, *t);
            to << "Vibrational:\t\t" << fix(evib) << '\t' << fix(cvib) << '\t'
               << fix(svib) << '\n';

            if (mol.has_torsions()) {
                double etor = chem::thermal_energy_tor(mol, *t) / datum::kilo;
                double ctor = chem::const_vol_heat_tor(mol, *t);
                double stor = chem::entropy_tor(mol, *t);
                to << "Torsional:\t\t" << fix(etor) << '\t' << fix(ctor) << '\t'
                   << fix(stor) << '\n';
            }

            line.width(36).fill('-');
            to << "\n\t\t\t"
               << "Q(" << *t << " K)\n"
               << line('-') << '\n';

            chem::Format<double> sci;
            sci.scientific().width(12).precision(6);

            double q = chem::qtot(mol, *t, *p, incl_sigma, "BOT");
            to << "Total (BOT):\t\t" << sci(q) << '\n';

            q = chem::qtot(mol, *t, *p, incl_sigma, "V=0");
            to << "Total (V=0):\t\t" << sci(q) << '\n';

            q = chem::qvib(mol, *t, "BOT");
            to << "Vibr. (BOT):\t\t" << sci(q) << '\n';

            q = chem::qvib(mol, *t, "V=0");
            to << "Vibr. (V=0):\t\t" << sci(q) << '\n';

            q = chem::qelec(mol, *t);
            to << "Electronic:\t\t" << sci(q) << '\n';

            q = chem::qtrans(mol, *t, *p);
            to << "Translational:\t\t" << sci(q) << '\n';

            q = chem::qrot(mol, *t, incl_sigma);
            to << "Rotational:\t\t" << sci(q) << '\n';

            if (mol.has_torsions()) {
                q = chem::qtor(mol, *t);
                to << "Torsional:\t\t" << sci(q) << '\n';
            }
            to << '\n';
        }
    }
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