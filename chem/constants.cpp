/**
   @file constants.cpp
   
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

#include <chem/constants.h>


namespace constants {

    const double m_e        = 2.7182818284590452354;
    const double m_log2e    = 1.4426950408889634074;
    const double m_log10e   = 0.43429448190325182765;
    const double m_ln2      = 0.69314718055994530942;
    const double m_ln10     = 2.30258509299404568402;
    const double m_pi       = 3.14159265358979323846;
    const double m_twopi    = 2.0 * m_pi;
    const double m_pi_2     = 1.57079632679489661923;
    const double m_pi_4     = 0.78539816339744830962;
    const double m_3pi_4    = 2.3561944901923448370e0;
    const double m_sqrtpi   = 1.77245385090551602792981;
    const double m_1_pi     = 0.31830988618379067154;
    const double m_2_pi     = 0.63661977236758134308;
    const double m_2_sqrtpi = 1.12837916709551257390;
    const double m_sqrt2    = 1.41421356237309504880;
    const double m_sqrt1_2  = 0.70710678118654752440;
    const double m_sqrt3    = 1.73205080756887719000;
    const double m_invln10  = 0.43429448190325182765;
    const double m_invln2   = 1.4426950408889633870e0;

    const double atomic_mass      = 1.66053904000e-27;
    const double avogadro         = 6.02214085700e+23;
    const double bohr_radius      = 0.529177210670;
    const double boltzmann        = 1.38064852000e-23;   
    const double conduct_quant    = 7.74809173100e-5;
    const double elec_const       = 8.85418781700e-12;
    const double electron_mass    = 9.10938356000e-31;
    const double electron_volt    = 1.60217662080e-19;
    const double elem_charge      = 1.60217662080e-19;
    const double faraday          = 9.64853328900e+4;
    const double fine_structure   = 7.29735256640e-3;
    const double gas_const        = 8.31445980000;        
    const double hartree          = 4.35974465000e-18;
    const double light_speed      = 2.99792458e+8;   
    const double mag_const        = 1.25663706140e-6;
    const double mag_flux_quant   = 2.06783383100e-15;
    const double mp_me_ratio      = 1.83615267389e+3;
    const double newtonian_grav   = 6.67408000000e-11;
    const double planck           = 6.62607004000e-34;  
    const double planck_bar       = 1.05457180000e-34;  
    const double proton_mass      = 1.67262189800e-27;
    const double rydberg          = 1.09737315685e+7;
    const double std_atm          = 1.01325000000e+5;        
    const double stefan_boltzmann = 5.67036700000e-8;

    const double yotta = 1.0e+24;
    const double zetta = 1.0e+21;
    const double exa   = 1.0e+18;
    const double peta  = 1.0e+15;
    const double tera  = 1.0e+12;
    const double giga  = 1.0e+9;
    const double mega  = 1.0e+6;
    const double kilo  = 1.0e+3;
    const double hecto = 1.0e+2;
    const double deca  = 10.0;
    const double one   = 1.0;
    const double deci  = 1.0e-1;
    const double centi = 1.0e-2;
    const double milli = 1.0e-3;
    const double micro = 1.0e-6;
    const double nano  = 1.0e-9;
    const double pico  = 1.0e-12;
    const double femto = 1.0e-15;
    const double atto  = 1.0e-18;
    const double zepto = 1.0e-21;
    const double yocto = 1.0e-24;

    const double cal2J  = 4.184;          
    const double icm2kJ = 1.19626564e-02; 
    const double icm2K  = 100.0 * planck * light_speed / boltzmann;
    const double J2icm  = 100.0 * planck * light_speed;
    const double au2cm  = bohr_radius * 1.0e-8;
    const double au2icm = hartree / (planck * light_speed * 100.0);
    const double au2s   = planck_bar / hartree;
    const double au2K   = hartree / boltzmann;
    const double au2kg  = electron_mass;

} // constants::

