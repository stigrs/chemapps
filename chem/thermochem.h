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
#include <cmath>
#include <gsl/gsl>

namespace chem {

// Translational:

// Calculate translational partitition function.
double qtrans(const Molecule& mol, double temp = 298.15, double pressure = 0.0);

// Calculate translational contribution to entropy.
double entropy_trans(const Molecule& mol,
                     double temp     = 298.15,
                     double pressure = datum::std_atm);

// Calculate translational contribution to internal thermal energy:
double thermal_energy_trans(double temp = 298.15);

// Calculate translational constant volume heat capacity.
double const_vol_heat_trans() { return 1.5 * datum::R; }

// Electronic:

// Calculate electronic partition function.
double qelec(const Molecule& mol, double temp = 298.15);

// Calculate electronic contribution to entropy.
double entropy_elec(const Molecule& mol, double temp = 298.15);

// Calculate electronic contribution to internal thermal energy.
double thermal_energy_elec() { return 0.0; }

// Calculate electronic constant volume heat capacity.
double const_vol_heat_elec() { return 0.0; }

// Rotational:

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

// Vibrational:

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

// Torsional:

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

inline double chem::entropy_elec(const Molecule& mol, double temp)
{
    double qe = chem::qelec(mol, temp);
    Ensures(qe > 0.0);
    return datum::R * std::log(qe);
}

inline double chem::thermal_energy_rot(const Molecule& mol, double temp)
{
    Expects(temp >= 0.0);
    return chem::const_vol_heat_rot(mol) * temp;
}

inline double chem::const_vol_heat_rot(const Molecule& mol)
{
    std::string rot_symm = mol.get_rot().symmetry();

    double cv_r = 0.0;
    if (rot_symm.find("atom") != std::string::npos) {
        cv_r = 0.0;
    }
    else {
        double factor = 1.5;
        if (rot_symm.find("linear") != std::string::npos) {
            factor = 1.0;
        }
        cv_r = factor * datum::R;
    }
    return cv_r;
}

#endif  // CHEM_THERMOCHEM_H
