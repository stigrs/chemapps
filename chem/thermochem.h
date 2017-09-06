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

#ifndef CHEM_THERMOCHEM_H
#define CHEM_THERMOCHEM_H

#include <chem/datum.h>
#include <chem/molecule.h>
#include <chem/utils.h>
#include <armadillo>
#include <cmath>
#include <gsl/gsl>
#include <iostream>
#include <limits>

namespace chem {

// Perform thermochemistry analysis.
void thermochemistry(const Molecule& mol,
                     const arma::vec& temp     = arma::vec{298.15},
                     const arma::vec& pressure = arma::vec{datum::std_atm},
                     bool incl_sigma           = true,
                     std::ostream& to          = std::cout);

//
// Translational:
//

// Calculate translational partitition function.
double qtrans(const Molecule& mol, double temp = 298.15, double pressure = 0.0);

// Calculate translational contribution to entropy.
double entropy_trans(const Molecule& mol,
                     double temp     = 298.15,
                     double pressure = datum::std_atm);

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
double entropy_rot(const Molecule& mol,
                   double temp     = 298.15,
                   bool incl_sigma = true);

// Calculate rotational contribution ot internal thermal energy.
double thermal_energy_rot(const Molecule& mol, double temp = 298.15);

// Calculate rotational contribution to constant volume heat capacity.
double const_vol_heat_rot(const Molecule& mol);

//
// Vibrational:
//

// Calculate vibrational partition function.
double qvib(const Molecule& mol,
            double temp                = 298.15,
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
            double temp               = 298.15,
            const std::string& scheme = "CT-Cw");

// Calculate partition function for a torsional mode using the CT-Cw scheme.
// Chuang, Y. Y.; Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.
double qctcw(const Molecule& mol, double temp = 298.15);

// Function for computing the derivative dln(Q)/dT for torsional modes.
double dlnqtor_dt(const Molecule& mol,
                  double temp               = 298.15,
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
            double temp                = 298.15,
            double pressure            = datum::std_atm,
            bool incl_sigma            = true,
            const std::string& zeroref = "BOT");

// Calculate total entropy.
double entropy(const Molecule& mol,
               double temp     = 298.15,
               double pressure = datum::std_atm,
               bool incl_sigma = true);

// Calculate thermal correction to the energy.
double thermal_energy(const Molecule& mol, double temp = 298.15);

// Calculate total constant volume heat capacity.
double const_vol_heat_capacity(const Molecule& mol, double temp = 298.15);

// Calculate thermal correction to enthalpy.
double enthalpy(const Molecule& mol, double temp = 298.15);

// Calculate thermal correction to Gibbs energy.
double gibbs_energy(const Molecule& mol,
                    double temp     = 298.15,
                    double pressure = datum::std_atm,
                    bool incl_sigma = true);

}  // namespace chem

inline double chem::qtrans(const Molecule& mol, double temp, double pressure)
{
    using namespace datum;

    Expects(temp >= 0.0);
    Expects(pressure >= 0.0);

    double mass = mol.tot_mass() * m_u;
    double vol  = 1.0;
    if (pressure > 0.0) {
        vol = k * temp / pressure;
    }
    return std::pow((2.0 * pi * mass * k * temp) / (h * h), 1.5) * vol;
}

inline double chem::entropy_trans(const Molecule& mol,
                                  double temp,
                                  double pressure)
{
    double qt = chem::qtrans(mol, temp, pressure);
    Ensures(qt > 0.0);
    return datum::R * (std::log(qt) + 2.5);
}

inline double chem::thermal_energy_trans(double temp)
{
    Expects(temp >= 0.0);
    return 1.5 * datum::R * temp;
}

inline double chem::const_vol_heat_trans() { return 1.5 * datum::R; }

inline double chem::entropy_elec(const Molecule& mol, double temp)
{
    double qe = chem::qelec(mol, temp);
    Ensures(qe > 0.0);
    return datum::R * std::log(qe);
}

inline double chem::thermal_energy_elec() { return 0.0; }

inline double chem::const_vol_heat_elec() { return 0.0; }

inline double chem::thermal_energy_rot(const Molecule& mol, double temp)
{
    Expects(temp >= 0.0);
    return chem::const_vol_heat_rot(mol) * temp;
}

inline double chem::const_vol_heat_rot(const Molecule& mol)
{
    std::string rot_symm = mol.get_rot().symmetry();

    if (rot_symm.find("atom") != std::string::npos) {
        return 0.0;
    }
    else {
        double factor = 1.5;
        if (rot_symm.find("linear") != std::string::npos) {
            factor = 1.0;
        }
        return factor * datum::R;
    }
}

inline double chem::qtor(const Molecule& mol,
                         double temp,
                         const std::string& scheme)
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

inline double chem::dlnqtor_dt(const Molecule& mol,
                               double temp,
                               const std::string& scheme)
{
    // The derivative is computed numerically using central difference.

    Expects(temp > 0.0);
    double h =  // Numerical recipes, Ch. 5.7
        std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0) * temp;
    double lnqa = std::log(chem::qtor(mol, temp + h, scheme));
    double lnqb = std::log(chem::qtor(mol, temp - h, scheme));
    return (lnqa - lnqb) / (2.0 * h);
}

inline double chem::entropy_tor(const Molecule& mol, double temp)
{
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
            return datum::R
                   * (std::log(chem::qtor(mol, temp))
                      + temp * chem::dlnqtor_dt(mol, temp));
        }
    }
}

inline double chem::thermal_energy_tor(const Molecule& mol, double temp)
{
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
            return datum::R * temp * temp * chem::dlnqtor_dt(mol, temp);
        }
    }
}

inline double chem::qtot(const Molecule& mol,
                         double temp,
                         double pressure,
                         bool incl_sigma,
                         const std::string& zeroref)
{
    return chem::qelec(mol, temp) * chem::qtrans(mol, temp, pressure)
           * chem::qrot(mol, temp, incl_sigma) * chem::qvib(mol, temp, zeroref)
           * chem::qtor(mol, temp);
}

inline double chem::entropy(const Molecule& mol,
                            double temp,
                            double pressure,
                            bool incl_sigma)
{
    return chem::entropy_elec(mol, temp)
           + chem::entropy_trans(mol, temp, pressure)
           + chem::entropy_rot(mol, temp, incl_sigma)
           + chem::entropy_vib(mol, temp) + chem::entropy_tor(mol, temp);
}

inline double chem::thermal_energy(const Molecule& mol, double temp)
{
    return chem::thermal_energy_elec() + chem::thermal_energy_trans(temp)
           + chem::thermal_energy_rot(mol, temp)
           + chem::thermal_energy_vib(mol, temp)
           + chem::thermal_energy_tor(mol, temp);
}

inline double chem::const_vol_heat_capacity(const Molecule& mol, double temp)
{
    return chem::const_vol_heat_elec() + chem::const_vol_heat_trans()
           + chem::const_vol_heat_rot(mol) + chem::const_vol_heat_vib(mol, temp)
           + chem::const_vol_heat_tor(mol, temp);
}

inline double chem::enthalpy(const Molecule& mol, double temp)
{
    Expects(temp >= 0.0);
    return chem::thermal_energy(mol, temp) + datum::R * temp;
}

inline double chem::gibbs_energy(const Molecule& mol,
                                 double temp,
                                 double pressure,
                                 bool incl_sigma)
{
    Expects(temp >= 0.0);
    return chem::enthalpy(mol, temp)
           - temp * chem::entropy(mol, temp, pressure, incl_sigma);
}

#endif  // CHEM_THERMOCHEM_H
