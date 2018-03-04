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

#ifndef CHEM_COLLISION_H
#define CHEM_COLLISION_H

#include <chem/mol_formula.h>
#include <srs/datum.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Error reporting:

struct Collision_error : std::runtime_error {
    Collision_error(std::string s) : std::runtime_error(s) {}
};

//
// Class providing methods for computing collision integrals and Lennard-Jones
// collision rates.
//
// Algorithms:
//   Forst, W. Unimolecular Reactions; Cambridge University Press, 2003.
//   Gilbert, R. G., J. Chem. Phys., 1984, vol. 80, pp. 5501-5509.
//   Lim, K. F.; Gilbert, R. G. J. Chem. Phys., 1990, vol. 92, pp. 1819-1830.
//   Troe, J. J. Chem. Phys., 1977, vol. 66, pp. 4758-4775.
//
class Collision {
public:
    Collision(std::istream& from, const std::string& key = "Collision");

    // Calculate reduced mass of system.
    double reduced_mass() const;

    // Average atom/atom mass of molecule (Lim and Gilbert, 1990).
    double average_mass() const;

    // Calculate Lennard-Jones well depth of system.
    double epsilon_complex() const;

    // Calculate Lennard-Jones collision diameter of system.
    double sigma_complex() const;

    // Local Lennard-Jones well depth of system (Lim and Gilbert, 1990).
    double epsilon_local() const;

    // Local Lennard-Jones collision diam. of system (Lim and Gilbert, 1990).
    double sigma_local() const;

    // Calculate Lennard-Jones collision rate (Forst, 2003).
    double lj_coll_rate() const;

    // Calculate reduced collision integral.
    double coll_omega22() const;

    // Collision time using eq. 32 of Lim and Gilbert (1990).
    double collision_time() const;

    // Biased random walk parameter s using eq. 30 in Lim and Gilbert (1990).
    double s_parameter() const;

    // Mean-squared energy transfer per collision (<E^2>).
    double mean_sqr_energy_transfer_coll() const;

    // Calculate collision energy transfer in highly excited molecules
    // using the biased random walk model B.
    void biased_random_walk(std::ostream& to = std::cout) const;

private:
    // Populate database with local Lennard-Jones collision diameter values.
    void set_sigma_local_values();

    // Populate database with local Lennard-Jones well depth values.
    void set_epsilon_local_values();

    // Closest interaction distance (d) (Lim and Gilbert, 1990).
    double dist_interact() const;

    // Impact parameter (b) (Lim and Gilbert, 1990).
    double impact_parameter() const;

    // Average translational energy (Lim and Gilbert, 1990).
    double energy_trans_avg() const;

    // A decay parameter (Lim and Gilbert, 1990).
    double a_decay_parameter() const;

    // Autocorrelation oscillation frequency (Lim and Gilbert, 1990).
    double c_autocorr_osc_freq() const;

    // Mean-squared rate of internal energy change (Lim and Gilbert, 1990).
    double mean_sqr_int_energy_change() const;

    // Find lightest mass of molecule.
    double mol_mass_lightest() const;

    enum Coll_omega22_t { troe, forst };  // collision integral equations

    Coll_omega22_t coll_integral;

    double mass_bath;     // mass of bath gas in amu
    double mass_mol;      // mass of molecule in amu
    double epsilon_bath;  // LJ well depth of bath gas in kelvin
    double epsilon_mol;   // LJ well depth of molecule in kelvin
    double sigma_bath;    // LJ collision diam. of bath gas in angstrom
    double sigma_mol;     // LJ collision diam. of molecule in angstrom
    double temperature;   // temperature in kelvin
    double vibr_high;     // highest vibrational frequency of molecule in cm-1

    std::vector<double> sigma_loc_val;     // local sigma values
    std::vector<double> epsilon_loc_val;   // local epsilon values
    std::vector<Mol_formula> mol_formula;  // molecular formula of collider
};

inline double Collision::reduced_mass() const
{
    return mass_bath * mass_mol / (mass_bath + mass_mol);
}

inline double Collision::epsilon_complex() const
{
    return std::sqrt(epsilon_bath * epsilon_mol);
}

inline double Collision::sigma_complex() const
{
    return 0.5 * (sigma_bath + sigma_mol);
}

inline double Collision::lj_coll_rate() const
{
    double sig = sigma_complex();
    double mu  = reduced_mass();

    return 4.5713e-12 * sig * sig * std::sqrt(temperature / mu)
           * coll_omega22();
}

inline double Collision::s_parameter() const
{
    double edot = mean_sqr_int_energy_change();
    double tc   = collision_time();
    double a    = a_decay_parameter();
    double c    = c_autocorr_osc_freq();

    return std::sqrt(edot * tc * 2.0 * a / (a * a + c * c));
}

inline double Collision::mean_sqr_energy_transfer_coll() const
{
    double s = s_parameter();
    return 2.0 * s * s;
}

//------------------------------------------------------------------------------

inline double Collision::dist_interact() const
{
    double d = sigma_complex() * std::sqrt(coll_omega22());
    return std::max(d, sigma_complex());
}

inline double Collision::impact_parameter() const
{
    return (2.0 / 3.0) * dist_interact();
}

inline double Collision::energy_trans_avg() const
{
    return 2.0e-3 * datum::k * temperature * datum::N_A / datum::icm_to_kJ;
}

inline double Collision::c_autocorr_osc_freq() const
{
    const double nu = vibr_high * datum::c_0 * 100.0;
    return 2.0 * datum::pi * nu;
}

#endif  // CHEM_COLLISION_H
