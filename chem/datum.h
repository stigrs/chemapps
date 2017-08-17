/**
   @file datum.h

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

#ifndef CHEM_DATUM_H
#define CHEM_DATUM_H

/**
   Namespace providing mathematical and physical constants, metric prefixes,
   and conversion factors.
*/
namespace datum {

// Mathematical constants:
extern const double e;
extern const double pi;
extern const double sqrt2;

// CODATA recommended 2014 values for physical constants:
// Source: http://physics.nist.gov/constants
extern const double m_u;      ///< Atomic mass constant (kg)
extern const double N_A;      ///< Avogadro constant (mol^-1)
extern const double a_0;      ///< Bohr radius (1.0E-10 m)
extern const double k;        ///< Boltzmann constant (J K^-1)
extern const double G_0;      ///< Conductance quantum (S)
extern const double eps_0;    ///< Electric constant (F m^-1)
extern const double eV;       ///< Electron volt (J)
extern const double m_e;      ///< Electron mass (kg)
extern const double ec;       ///< Elementary charge (C)
extern const double F;        ///< Faraday constant (C mol^-1)
extern const double alpha;    ///< Fine-structure constant
extern const double R;        ///< Gas constant (J mol^-1 K^-1)
extern const double E_h;      ///< Hartree energy (J)
extern const double c_0;      ///< Speed of light (m s^-1)
extern const double mu_0;     ///< Magnetic constant (N A^-2)
extern const double phi_0;    ///< Magnetic flux quantum (Wb)
extern const double m_p_m_e;  ///< Proton-electron mass ratio
extern const double G;        ///< Newtonian grav. (m^3/(kg s^2))
extern const double h;        ///< Planck constant (J s)
extern const double h_bar;    ///< Planck-bar constant (J s)
extern const double m_p;      ///< Proton mass (kg)
extern const double R_inf;    ///< Rydberg constant (m^-1)
extern const double std_atm;  ///< Std. atmospheric pressure (Pa)
extern const double sigma;    ///< Stefan-Boltzmann (W m^-2 K^-4)

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
extern const double cal_to_J;    ///< Conversion from cal to J
extern const double icm_to_kJ;   ///< Conversion from cm^-1 to kJ/mol
extern const double icm_to_K;    ///< Conversion from cm^-1 to K
extern const double J_to_icm;    ///< Conversion from J to cm^-1
extern const double au_to_cm;    ///< Conversion from au to cm
extern const double au_to_icm;   ///< Conversion from au to cm^-1
extern const double au_to_s;     ///< Conversion from au to s
extern const double au_to_K;     ///< Conversion from au to K
extern const double au_to_kg;    ///< Conversion from au to kg
extern const double au_to_kgm2;  ///< Conversion from amu bohr^2 to kg m^2
}  // datum::

#endif /* CHEM_DATUM_H */
