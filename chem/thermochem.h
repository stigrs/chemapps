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
    double qtr = chem::qtrans(mol, temp, pressure);
    Ensures(qtr > 0.0);
    return datum::R * (std::log(qtr) + 2.5);
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
#endif  // CHEM_THERMOCHEM_H