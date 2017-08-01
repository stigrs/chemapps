/**
   @file constants.h
   
   This file is part of ChemApps - A C++ Chemistry Toolkit
   
   Copyright (C) 2016-2017  Stig Rune Sellevag
   
   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CHEM_CONSTANTS_H
#define CHEM_CONSTANTS_H


/**
   Namespace providing mathematical and physical constants, metric prefixes,
   and conversion factors.
*/
namespace constants {

    // Mathematical constants:
    extern const double m_e;
    extern const double m_log2e;    ///< log2(e)
    extern const double m_log10e;   ///< log10(e)
    extern const double m_ln2;      ///< ln(2)
    extern const double m_ln10;     ///< ln(10)
    extern const double m_pi;
    extern const double m_twopi;    ///< 2 pi
    extern const double m_pi_2;     ///< pi / 2
    extern const double m_pi_4;     ///< pi / 4
    extern const double m_3pi_4;    ///< 3pi / 4
    extern const double m_sqrtpi;   ///< sqrt(pi)
    extern const double m_1_pi;     ///< 1 / pi
    extern const double m_2_pi;     ///< 2 / pi
    extern const double m_2_sqrtpi; ///< 2 / sqrt(pi)
    extern const double m_sqrt2;    ///< sqrt(2)
    extern const double m_sqrt1_2;  ///< 1 / sqrt(2)
    extern const double m_sqrt3;    ///< sqrt(3)
    extern const double m_invln10;  ///< 1 / ln(10)
    extern const double m_invln2;   ///< 1 / ln(2)

    // CODATA recommended 2014 values for physical constants:
    // Source: http://physics.nist.gov/constants
    extern const double atomic_mass;      ///< Atomic mass constant (kg)
    extern const double avogadro;         ///< Avogadro constant (mol^-1)
    extern const double bohr_radius;      ///< Bohr radius (1.0E-10 m)
    extern const double boltzmann;        ///< Boltzmann constant (J K^-1)
    extern const double conduct_quant;    ///< Conductance quantum (S)
    extern const double elec_const;       ///< Electric constant (F m^-1)
    extern const double electron_volt;    ///< Electron volt (J)
    extern const double electron_mass;    ///< Electron mass (kg)
    extern const double elem_charge;      ///< Elementary charge (C)
    extern const double faraday;          ///< Faraday constant (C mol^-1)
    extern const double fine_structure;   ///< Fine-structure constant
    extern const double gas_const;        ///< Gas constant (J mol^-1 K^-1)
    extern const double hartree;          ///< Hartree energy (J)
    extern const double light_speed;      ///< Speed of light (m s^-1)
    extern const double mag_const;        ///< Magnetic constant (N A^-2)
    extern const double mag_flux_quant;   ///< Magnetic flux quantum (Wb)
    extern const double mp_me_ratio;      ///< Proton-electron mass ratio
    extern const double newtonian_grav;   ///< Newtonian grav. (m^3/(kg s^2))
    extern const double planck;           ///< Planck constant (J s)
    extern const double planck_bar;       ///< Planck-bar constant (J s)
    extern const double proton_mass;      ///< Proton mass (kg)
    extern const double rydberg;          ///< Rydberg constant (m^-1)
    extern const double std_atm;          ///< Std. atmospheric pressure (Pa)
    extern const double stefan_boltzmann; ///< Stefan-Boltzmann (W m^-2 K^-4)

    // Metric prefixes:
    extern const double yotta;
    extern const double zetta;
    extern const double exa;
    extern const double peta;
    extern const double tera;
    extern const double giga;
    extern const double mega;
    extern const double kilo;
    extern const double hecto;
    extern const double deca;
    extern const double one;
    extern const double deci;
    extern const double centi;
    extern const double milli;
    extern const double micro;
    extern const double nano;
    extern const double pico;
    extern const double femto;
    extern const double atto;
    extern const double zepto;
    extern const double yocto;

    // Conversion factors:
    extern const double cal2J;  ///< Conversion factor from cal to J
    extern const double icm2kJ; ///< Conversion factor from cm^-1 to kJ/mol
    extern const double icm2K;  ///< Conversion factor from cm^-1 to K
    extern const double J2icm;  ///< Conversion factor from J to cm^-1
    extern const double au2cm;  ///< Conversion factor from au to cm
    extern const double au2icm; ///< Conversion factor from au to cm^-1
    extern const double au2s;   ///< Conversion factor from au to s
    extern const double au2K;   ///< Conversion factor from au to K
    extern const double au2kg;  ///< Conversion factor from au to kg
} // constants::

#endif /* CHEM_CONSTANTS_H */

