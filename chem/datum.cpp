//////////////////////////////////////////////////////////////////////////////
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

namespace datum {

const double e     = 2.7182818284590452354;
const double pi    = 3.14159265358979323846;
const double sqrt2 = 1.41421356237309504880;

const double m_u     = 1.66053904000e-27;
const double N_A     = 6.02214085700e+23;
const double a_0     = 0.529177210670;
const double k       = 1.38064852000e-23;
const double G_0     = 7.74809173100e-5;
const double eps_0   = 8.85418781700e-12;
const double m_e     = 9.10938356000e-31;
const double eV      = 1.60217662080e-19;
const double ec      = 1.60217662080e-19;
const double F       = 9.64853328900e+4;
const double alpha   = 7.29735256640e-3;
const double R       = 8.31445980000;
const double E_h     = 4.35974465000e-18;
const double c_0     = 2.99792458e+8;
const double mu_0    = 1.25663706140e-6;
const double phi_0   = 2.06783383100e-15;
const double m_p_m_e = 1.83615267389e+3;
const double G       = 6.67408000000e-11;
const double h       = 6.62607004000e-34;
const double h_bar   = 1.05457180000e-34;
const double m_p     = 1.67262189800e-27;
const double R_inf   = 1.09737315685e+7;
const double std_atm = 1.01325000000e+5;
const double sigma   = 5.67036700000e-8;

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

const double cal_to_J   = 4.184;
const double icm_to_kJ  = 1.19626564e-02;
const double icm_to_K   = 100.0 * h * c_0 / k;
const double J_to_icm   = 100.0 * h * c_0;
const double au_to_cm   = a_0 * 1.0e-8;
const double au_to_icm  = E_h / (h * c_0 * 100.0);
const double au_to_s    = h_bar / E_h;
const double au_to_K    = E_h / k;
const double au_to_kg   = m_e;
const double au_to_kgm2 = m_u * a_0 * a_0 * 1.0e-20;

}  // namespace datum
