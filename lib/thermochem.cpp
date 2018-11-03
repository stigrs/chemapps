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

#include <chem/thermochem.h>
#include <numlib/math.h>
#include <string>

void Chem::thermochemistry(const Chem::Molecule& mol,
                           const Numlib::Vec<double>& temp,
                           const Numlib::Vec<double>& pressure,
                           bool incl_sigma,
                           std::ostream& to)
{
    using namespace Numlib;

    Chem::Molecule tmp(mol);

    Stdutils::Format<char> line;
    if (!tmp.info().empty()) {
        line.width(20 + mol.info().size()).fill('=');
    }
    else {
        line.width(16).fill('=');
    }
    Stdutils::Format<double> fix;
    fix.fixed().precision(6);

    double e0 = tmp.elec_energy();
    if (!tmp.info().empty()) {
        to << "\nThermochemistry of " << tmp.info() << ":\n"
           << line('=') << '\n';
    }
    else {
        to << "\nThermochemistry:\n" << line('=') << '\n';
    }
    to << "Electronic energy: " << fix(e0) << " Hartree\n";

    tmp.rot_analysis(to);

    if (tmp.tot_tor_minima() > 0) {
        tmp.tor_analysis(to);
    }
    tmp.vib_analysis(to);

    double zpe = tmp.zero_point_energy() / Constants::au_to_icm;
    const double factor = 1.0 / (Constants::E_h * Constants::N_A);

    for (const auto& p : pressure) {
        for (const auto& t : temp) {
            to << "Temperature: " << t << " K. "
               << "Pressure: " << p << " bar\n";

            fix.fixed().width(12).precision(6);
            to << "Zero-point correction:\t\t\t\t" << fix(zpe) << " Hartree\n";

            double ecorr = Chem::thermal_energy(tmp, t) * factor;
            to << "Thermal correction to energy:\t\t\t" << fix(ecorr) << '\n';

            double hcorr = Chem::enthalpy(tmp, t) * factor;
            to << "Thermal correction to enthalpy:\t\t\t" << fix(hcorr) << '\n';

            double gcorr = Chem::gibbs_energy(tmp, t, p, incl_sigma) * factor;
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
            etot = ecorr * Constants::au_to_icm * Constants::icm_to_kJ;
            double cv = Chem::const_vol_heat_capacity(tmp, t);
            double stot = Chem::entropy(tmp, t, p, incl_sigma);
            to << "Total:\t\t\t" << fix(etot) << '\t' << fix(cv) << '\t'
               << fix(stot) << '\n';

            double eelec = Chem::thermal_energy_elec();
            double celec = Chem::const_vol_heat_elec();
            double selec = Chem::entropy_elec(tmp, t);
            to << "Electronic:\t\t" << fix(eelec) << '\t' << fix(celec) << '\t'
               << fix(selec) << '\n';

            double etrans = Chem::thermal_energy_trans(t) / Constants::kilo;
            double ctrans = Chem::const_vol_heat_trans();
            double strans = Chem::entropy_trans(tmp, t, p);
            to << "Translational:\t\t" << fix(etrans) << '\t' << fix(ctrans)
               << '\t' << fix(strans) << '\n';

            double erot = Chem::thermal_energy_rot(tmp, t) / Constants::kilo;
            double crot = Chem::const_vol_heat_rot(tmp);
            double srot = Chem::entropy_rot(tmp, t, incl_sigma);
            to << "Rotational:\t\t" << fix(erot) << '\t' << fix(crot) << '\t'
               << fix(srot) << '\n';

            double evib = Chem::thermal_energy_vib(tmp, t) / Constants::kilo;
            double cvib = Chem::const_vol_heat_vib(tmp, t);
            double svib = Chem::entropy_vib(tmp, t);
            to << "Vibrational:\t\t" << fix(evib) << '\t' << fix(cvib) << '\t'
               << fix(svib) << '\n';

            if (tmp.tot_tor_minima() > 0) {
                double etor =
                    Chem::thermal_energy_tor(tmp, t) / Constants::kilo;
                double ctor = Chem::const_vol_heat_tor(tmp, t);
                double stor = Chem::entropy_tor(tmp, t);
                to << "Torsional:\t\t" << fix(etor) << '\t' << fix(ctor) << '\t'
                   << fix(stor) << '\n';
            }

            line.width(36).fill('-');
            to << "\n\t\t\t"
               << "Q(" << t << " K)\n"
               << line('-') << '\n';

            Stdutils::Format<double> sci;
            sci.scientific().width(12).precision(6);

            double q = Chem::qtot(tmp, t, p, incl_sigma, "BOT");
            to << "Total (BOT):\t\t" << sci(q) << '\n';

            q = Chem::qtot(tmp, t, p, incl_sigma, "V=0");
            to << "Total (V=0):\t\t" << sci(q) << '\n';

            q = Chem::qvib(tmp, t, "BOT");
            to << "Vibr. (BOT):\t\t" << sci(q) << '\n';

            q = Chem::qvib(tmp, t, "V=0");
            to << "Vibr. (V=0):\t\t" << sci(q) << '\n';

            q = Chem::qelec(tmp, t);
            to << "Electronic:\t\t" << sci(q) << '\n';

            q = Chem::qtrans(tmp, t, p);
            to << "Translational:\t\t" << sci(q) << '\n';

            q = Chem::qrot(tmp, t, incl_sigma);
            to << "Rotational:\t\t" << sci(q) << '\n';

            if (tmp.tot_tor_minima() > 0) {
                q = Chem::qtor(tmp, t);
                to << "Torsional:\t\t" << sci(q) << '\n';
            }
            to << '\n';
        }
    }
}

double Chem::qelec(const Chem::Molecule& mol, double temp)
{
    using namespace Numlib::Constants;

    auto so_degen = mol.spin_orbit_degen();
    auto so_energy = mol.spin_orbit_energy();

    double qe = 0.0;
    for (Index i = 0; i < so_degen.size(); ++i) {
        qe += so_degen(i) * std::exp(-so_energy(i) * icm_to_K / temp);
    }
    return qe;
}

double Chem::qrot(const Chem::Molecule& mol, double temp, bool incl_sigma)
{
    using namespace Numlib::Constants;

    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");

    Chem::Molecule tmp(mol);

    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 1.0;
    }
    else if (rot_symm.find("linear") != std::string::npos) {
        auto rotc = GHz_to_K * tmp.rot_constants();
        double rsig = 1.0;
        if (incl_sigma) {
            rsig /= tmp.rot_sigma();
        }
        res = rsig * temp / rotc(0);
    }
    else { // nonlinear molecule
        auto rotc = GHz_to_K * tmp.rot_constants();
        double b = Numlib::prod(rotc);
        double rsig = std::sqrt(pi);
        if (incl_sigma) {
            rsig /= tmp.rot_sigma();
        }
        res = rsig * std::pow(temp, 1.5) / std::sqrt(b);
    }
    return res;
}

double
Chem::entropy_rot(const Chem::Molecule& mol, double temp, bool incl_sigma)
{
    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 1.0;
    }
    else {
        double factor = 1.5;
        if (rot_symm.find("linear") != std::string::npos) {
            factor = 1.0;
        }
        double qr = Chem::qrot(tmp, temp, incl_sigma);
        Assert::dynamic<Assert::level(2)>(qr > 0.0);
        res = Numlib::Constants::R * (std::log(qr) + factor);
    }
    return res;
}

double
Chem::qvib(const Chem::Molecule& mol, double temp, const std::string& zeroref)
{
    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 1.0;
    }
    else {
        Assert::dynamic<Assert::level(2)>(temp > 0.0);
        auto w = Numlib::Constants::icm_to_K * tmp.frequencies();
        double qv = 1.0;
        if (zeroref == "V=0") { // zero at first vibrational level
            for (auto wi : w) {
                if (wi < 0.0) { // ignore imaginary frequencies
                    continue;
                }
                else {
                    qv /= (1.0 - std::exp(-wi / temp));
                }
            }
        }
        else if (zeroref == "BOT") { // zero at the bottom of the well
            for (auto wi : w) {
                if (wi < 0.0) {
                    continue;
                }
                else {
                    qv *= std::exp(-wi / (2.0 * temp)) /
                          (1.0 - std::exp(-wi / temp));
                }
            }
        }
        res = qv;
    }
    return res;
}

double Chem::entropy_vib(const Chem::Molecule& mol, double temp)
{
    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        Assert::dynamic<Assert::level(2)>(temp > 0.0);
        auto w = Numlib::Constants::icm_to_K * tmp.frequencies();
        double sv = 0.0;
        for (auto wi : w) {
            if (wi < 0.0) { // ignore imaginary frequencies
                continue;
            }
            else {
                double wt = wi / temp;
                sv += wt / (std::exp(wt) - 1.0) - std::log(1.0 - std::exp(-wt));
            }
        }
        sv *= Numlib::Constants::R;
        res = sv;
    }
    return res;
}

double Chem::thermal_energy_vib(const Chem::Molecule& mol, double temp)
{
    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        Assert::dynamic<Assert::level(2)>(temp > 0.0);
        auto w = Numlib::Constants::icm_to_K * tmp.frequencies();
        double ev = 0.0;
        for (auto wi : w) {
            if (wi < 0.0) {
                continue;
            }
            else {
                ev += wi * (0.5 + 1.0 / (std::exp(wi / temp) - 1.0));
            }
        }
        ev *= Numlib::Constants::R;
        res = ev;
    }
    return res;
}

double Chem::const_vol_heat_vib(const Chem::Molecule& mol, double temp)
{
    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        Assert::dynamic<Assert::level(2)>(temp > 0.0);
        auto w = Numlib::Constants::icm_to_K * tmp.frequencies();
        double cv_v = 0.0;
        for (auto wi : w) {
            if (wi < 0.0) { // ignore imaginary frequencies
                continue;
            }
            else {
                double wt = wi / temp;
                cv_v +=
                    wt * wt * std::exp(wt) / std::pow(std::exp(wt) - 1.0, 2.0);
            }
        }
        cv_v *= Numlib::Constants::R;
        res = cv_v;
    }
    return res;
}

double Chem::qctcw(const Chem::Molecule& mol, double temp)
{
    using namespace Numlib::Constants;

    double qtor = 1.0;

    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    if (rot_symm.find("atom") != std::string::npos) {
        qtor = 1.0;
    }
    else {
        if (tmp.tot_tor_minima() > 0) {
            Assert::dynamic<Assert::level(2)>(temp > 0.0);
            // Calculate free rotor partition function:
            double imom = tmp.tor_eff_moment();
            imom *= au_to_kgm2;
            double sig = tmp.tor_symmetry_number();
            double qfr = std::sqrt(2.0 * pi * imom * k * temp) / (h_bar * sig);

            // Calculate partition function for harmonic oscillator and
            // intermediate case:
            double qho = 0.0;
            double qin = 0.0;
            auto pot = tmp.tor_pot_coeff();
            auto freq = tmp.tor_frequencies();
            Assert::dynamic<Assert::level(2)>(pot.size() == freq.size());
            for (Index i = 0; i < pot.size(); ++i) {
                double ui = pot(i) * icm_to_K;
                double wi = freq(i) * icm_to_K;
                qho += std::exp(-(ui + 0.5 * wi) / temp) /
                       (1.0 - std::exp(-wi / temp));
                qin += std::exp(-ui / temp) / (wi / temp);
            }
            qtor = qho * std::tanh(qfr / qin); // eq. 11 in C&T (2000).
        }
    }
    return qtor;
}

double Chem::const_vol_heat_tor(const Chem::Molecule& mol, double temp)
{
    // The constant volume heat capacity is calculated as the numerical
    // derivative of the thermal torsional energy with respect to
    // temperature (dEtor/dT) at constant N and V.

    Chem::Molecule tmp(mol);
    std::string rot_symm = tmp.rot_symmetry();
    double res = 0.0;

    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        if (rot_symm.find("linear") != std::string::npos) {
            res = 0.0; // a linear molecule cannot have torsional modes
        }
        else {
            Assert::dynamic<Assert::level(2)>(temp > 0.0);
            double h = temp * std::pow(std::numeric_limits<double>::epsilon(),
                                       1.0 / 3.0);
            double cva = Chem::thermal_energy_tor(tmp, temp + h);
            double cvb = Chem::thermal_energy_tor(tmp, temp - h);
            res = (cva - cvb) / (2.0 * h);
        }
    }
    return res;
}

