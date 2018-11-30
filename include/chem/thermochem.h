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

#ifndef CHEM_THERMOCHEM_H
#define CHEM_THERMOCHEM_H

#include <chem/molecule.h>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <cmath>
#include <iostream>
#include <limits>

namespace Chem {

// Perform thermochemistry analysis.
void thermochemistry(
    const Molecule& mol,
    const Numlib::Vec<double>& temp = Numlib::Vec<double>{298.15},
    const Numlib::Vec<double>& pressure =
        Numlib::Vec<double>{Numlib::Constants::std_atm},
    bool incl_sigma = true,
    std::ostream& to = std::cout);

//
// Translational:
//

// Calculate translational partitition function.
double qtrans(const Molecule& mol, double temp = 298.15, double pressure = 0.0);
double qtrans(double mass, double temp = 298.15);

// Calculate translational contribution to entropy.
double entropy_trans(const Molecule& mol,
                     double temp = 298.15,
                     double pressure = Numlib::Constants::std_atm);

// Calculate translational contribution to internal thermal energy:
double thermal_energy_trans(double temp = 298.15);

// Calculate translational constant volume heat capacity.
double const_vol_heat_trans();

//
// Electronic:
//

// Calculate electronic partition function.
double qelec(const Molecule& mol, double temp = 298.15);

// Calculate electronic contribution to entropy.
double entropy_elec(const Molecule& mol, double temp = 298.15);

// Calculate electronic contribution to internal thermal energy.
double thermal_energy_elec();

// Calculate electronic constant volume heat capacity.
double const_vol_heat_elec();

//
// Rotational:
//

// Calculate rotational partition function.
double qrot(const Molecule& mol, double temp = 298.15, bool incl_sigma = true);

// Calculate rotational contribution to entropy.
double
entropy_rot(const Molecule& mol, double temp = 298.15, bool incl_sigma = true);

// Calculate rotational contribution ot internal thermal energy.
double thermal_energy_rot(const Molecule& mol, double temp = 298.15);

// Calculate rotational contribution to constant volume heat capacity.
double const_vol_heat_rot(const Molecule& mol);

//
// Vibrational:
//

// Calculate vibrational partition function.
double qvib(const Molecule& mol,
            double temp = 298.15,
            const std::string& zeroref = "BOT");

// Calculate vibrational contribution to entropy.
double entropy_vib(const Molecule& mol, double temp = 298.15);

// Calculate vibrational contribution to internal thermal energy.
double thermal_energy_vib(const Molecule& mol, double temp = 298.15);

// Calculate vibrational constant volume heat capacity.
double const_vol_heat_vib(const Molecule& mol, double temp = 298.15);

//
// Torsional:
//

// Calculate partition function for a molecular torsional mode.
double qtor(const Molecule& mol,
            double temp = 298.15,
            const std::string& scheme = "CT-Cw");

// Calculate partition function for a torsional mode using the CT-Cw scheme.
// Chuang, Y. Y.; Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.
double qctcw(const Molecule& mol, double temp = 298.15);

// Function for computing the derivative dln(Q)/dT for torsional modes.
double dlnqtor_dt(const Molecule& mol,
                  double temp = 298.15,
                  const std::string& scheme = "CT-Cw");

// Calculate torsional contribution to entropy.
double entropy_tor(const Molecule& mol, double temp = 298.15);

// Calculate torsional contribution to internal thermal energy.
double thermal_energy_tor(const Molecule& mol, double temp = 298.15);

// Calculate torsional constant volume heat capacity.
double const_vol_heat_tor(const Molecule& mol, double temp = 298.15);

//
// Total contributions:
//

// Calculate total molecular partition function.
double qtot(const Molecule& mol,
            double temp = 298.15,
            double pressure = Numlib::Constants::std_atm,
            bool incl_sigma = true,
            const std::string& zeroref = "BOT");

// Calculate total entropy.
double entropy(const Molecule& mol,
               double temp = 298.15,
               double pressure = Numlib::Constants::std_atm,
               bool incl_sigma = true);

// Calculate thermal correction to the energy.
double thermal_energy(const Molecule& mol, double temp = 298.15);

// Calculate total constant volume heat capacity.
double const_vol_heat_capacity(const Molecule& mol, double temp = 298.15);

// Calculate thermal correction to enthalpy.
double enthalpy(const Molecule& mol, double temp = 298.15);

// Calculate thermal correction to Gibbs energy.
double gibbs_energy(const Molecule& mol,
                    double temp = 298.15,
                    double pressure = Numlib::Constants::std_atm,
                    bool incl_sigma = true);

inline double qtrans(const Molecule& mol, double temp, double pressure)
{
    using namespace Numlib::Constants;

    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    Assert::dynamic<Assert::level(2)>(pressure >= 0.0, "bad pressure");

    double mass = mol.tot_mass() * m_u;
    double vol = 1.0;
    if (pressure > 0.0) {
        vol = k * temp / pressure;
    }
    return std::pow((2.0 * pi * mass * k * temp) / (h * h), 1.5) * vol;
}

inline double qtrans(double mass, double temp)
{
    using namespace Numlib::Constants;
    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    return std::pow(2.0 * pi * mass * m_u * k * temp, 1.5) / std::pow(h, 3.0);
}

inline double entropy_trans(const Molecule& mol, double temp, double pressure)
{
    double qt = qtrans(mol, temp, pressure);
    Assert::dynamic<Assert::level(2)>(qt > 0.0);
    return Numlib::Constants::R * (std::log(qt) + 2.5);
}

inline double thermal_energy_trans(double temp)
{
    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    return 1.5 * Numlib::Constants::R * temp;
}

inline double const_vol_heat_trans() { return 1.5 * Numlib::Constants::R; }

inline double entropy_elec(const Molecule& mol, double temp)
{
    double qe = qelec(mol, temp);
    Assert::dynamic<Assert::level(2)>(qe > 0.0);
    return Numlib::Constants::R * std::log(qe);
}

inline double thermal_energy_elec() { return 0.0; }

inline double const_vol_heat_elec() { return 0.0; }

inline double thermal_energy_rot(const Molecule& mol, double temp)
{
    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    return const_vol_heat_rot(mol) * temp;
}

inline double const_vol_heat_rot(const Molecule& mol)
{
    std::string rot_symm = mol.rot().symmetry();

    double res = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        double factor = 1.5;
        if (rot_symm.find("linear") != std::string::npos) {
            factor = 1.0;
        }
        res = factor * Numlib::Constants::R;
    }
    return res;
}

inline double qtor(const Molecule& mol, double temp, const std::string& scheme)
{
    // This is a wrapper function to allow future implementations of other
    // schemes for computing the partition function.

    if (scheme == "CT-Cw") {
        return qctcw(mol, temp);
    }
    else {
        return qctcw(mol, temp);
    }
}

inline double
dlnqtor_dt(const Molecule& mol, double temp, const std::string& scheme)
{
    // The derivative is computed numerically using central difference.

    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    double h = // Numerical recipes, Ch. 5.7
        std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0) * temp;
    double lnqa = std::log(qtor(mol, temp + h, scheme));
    double lnqb = std::log(qtor(mol, temp - h, scheme));
    return (lnqa - lnqb) / (2.0 * h);
}

inline double entropy_tor(const Molecule& mol, double temp)
{
    using namespace Numlib::Constants;

    std::string rot_symm = mol.rot().symmetry();

    double res = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        if (rot_symm.find("linear") != std::string::npos) {
            res = 0.0; // a linear molecule cannot have torsional modes
        }
        else {
            Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
            res =
                R * (std::log(qtor(mol, temp)) + temp * dlnqtor_dt(mol, temp));
        }
    }
    return res;
}

inline double thermal_energy_tor(const Molecule& mol, double temp)
{
    using namespace Numlib::Constants;

    std::string rot_symm = mol.rot().symmetry();

    double res = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        res = 0.0;
    }
    else {
        if (rot_symm.find("linear") != std::string::npos) {
            res = 0.0; // a linear molecule cannot have torsional modes
        }
        else {
            Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
            res = R * temp * temp * dlnqtor_dt(mol, temp);
        }
    }
    return res;
}

inline double qtot(const Molecule& mol,
                   double temp,
                   double pressure,
                   bool incl_sigma,
                   const std::string& zeroref)
{
    return qelec(mol, temp) * qtrans(mol, temp, pressure) *
           qrot(mol, temp, incl_sigma) * qvib(mol, temp, zeroref) *
           qtor(mol, temp);
}

inline double
entropy(const Molecule& mol, double temp, double pressure, bool incl_sigma)
{
    return entropy_elec(mol, temp) + entropy_trans(mol, temp, pressure) +
           entropy_rot(mol, temp, incl_sigma) + entropy_vib(mol, temp) +
           entropy_tor(mol, temp);
}

inline double thermal_energy(const Molecule& mol, double temp)
{
    return thermal_energy_elec() + thermal_energy_trans(temp) +
           thermal_energy_rot(mol, temp) + thermal_energy_vib(mol, temp) +
           thermal_energy_tor(mol, temp);
}

inline double const_vol_heat_capacity(const Molecule& mol, double temp)
{
    return const_vol_heat_elec() + const_vol_heat_trans() +
           const_vol_heat_rot(mol) + const_vol_heat_vib(mol, temp) +
           const_vol_heat_tor(mol, temp);
}

inline double enthalpy(const Molecule& mol, double temp)
{
    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    return thermal_energy(mol, temp) + Numlib::Constants::R * temp;
}

inline double
gibbs_energy(const Molecule& mol, double temp, double pressure, bool incl_sigma)
{
    Assert::dynamic<Assert::level(2)>(temp >= 0.0, "bad temperature");
    return enthalpy(mol, temp) -
           temp * entropy(mol, temp, pressure, incl_sigma);
}

} // namespace Chem

#endif // CHEM_THERMOCHEM_H

