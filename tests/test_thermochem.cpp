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

#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <chem/thermodata.h>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_thermochem")
{
    SECTION("CH3OH")
    {
        using namespace Numlib::Constants;

        std::ifstream from;
        Stdutils::fopen(from, "test_ch3oh.inp");
        Chem::Molecule mol(from);

        const double qtr_ans = 0.712383e+7;
        double qtr = Chem::qtrans(mol, 298.15, std_atm);
        CHECK(std::abs(qtr - qtr_ans) / qtr_ans < 5.0e-6);

        const double htr_ans = 0.889 * cal_to_J * 1000.0;
        double htr = Chem::thermal_energy_trans();
        CHECK(std::abs(htr - htr_ans) / htr_ans < 5.0e-4);

        const double str_ans = 36.324 * cal_to_J;
        double str = Chem::entropy_trans(mol);
        CHECK(std::abs(str - str_ans) / str_ans < 1.0e-8);

        const double cv_tr_ans = 2.981 * cal_to_J;
        double cv_tr = Chem::const_vol_heat_trans();
        CHECK(std::abs(cv_tr - cv_tr_ans) / cv_tr_ans < 1.0e-4);

        const double qe_ans = 0.100000e+1;
        double qe = Chem::qelec(mol);
        CHECK(std::abs(qe - qe_ans) / qe_ans < 1.0e-8);

        const double he_ans = 0.0;
        double he = Chem::thermal_energy_elec();
        CHECK(std::abs(he - he_ans) < 1.0e-9);

        const double cv_e_ans = 0.0;
        double cv_e = Chem::const_vol_heat_elec();
        CHECK(std::abs(cv_e - cv_e_ans) < 1.0e-9);

        const double qr_ans = 0.314134e+4;
        double qr = Chem::qrot(mol);
        CHECK(std::abs(qr - qr_ans) / qr_ans < 1.0e-5);

        const double sr_ans = 18.983 * cal_to_J;
        double sr = Chem::entropy_rot(mol);
        CHECK(std::abs(sr - sr_ans) / sr_ans < 3.0e-5);

        const double er_ans = 0.889 * cal_to_J * 1000.0;
        double er = Chem::thermal_energy_rot(mol);
        CHECK(std::abs(er - er_ans) / er_ans < 5.0e-4);

        const double cv_r_ans = 2.981 * cal_to_J;
        double cv_r = Chem::const_vol_heat_rot(mol);
        CHECK(std::abs(cv_r - cv_r_ans) / cv_r_ans < 1.0e-4);

        const double qv_ans = 0.148820e-23;
        double qv = Chem::qvib(mol);
        CHECK(std::abs(qv - qv_ans) / qv_ans < 1.0e-4);

        const double sv_ans = 1.441 * cal_to_J;
        double sv = Chem::entropy_vib(mol);
        CHECK(std::abs(sv - sv_ans) / sv_ans < 2.0e-4);

        const double ev_ans = 32.936 * cal_to_J * 1000.0;
        double ev = Chem::thermal_energy_vib(mol);
        CHECK(std::abs(ev - ev_ans) / ev_ans < 5.0e-5);

        const double cv_v_ans = 2.696 * cal_to_J;
        double cv_v = Chem::const_vol_heat_vib(mol);
        CHECK(std::abs(cv_v - cv_v_ans) / cv_v_ans < 2.0e-4);

        const double qtot_ans = 0.333036e-13;
        double qtot = Chem::qtot(mol);
        CHECK(std::abs(qtot - qtot_ans) / qtot_ans < 1.0e-4);

        const double stot_ans = 56.748 * cal_to_J;
        double stot = Chem::entropy(mol);
        CHECK(std::abs(stot - stot_ans) / stot_ans < 1.0e-5);

        const double etot_ans = 34.714 * cal_to_J * 1000.0;
        double etot = Chem::thermal_energy(mol);
        CHECK(std::abs(etot - etot_ans) / etot_ans < 2.0e-5);

        const double cv_tot_ans = 8.657 * cal_to_J;
        double cv_tot = Chem::const_vol_heat_capacity(mol);
        CHECK(std::abs(cv_tot - cv_tot_ans) / cv_tot_ans < 3.0e-5);

        const double hcorr_ans = 35.3061979 * cal_to_J * 1000.0;
        double hcorr = Chem::enthalpy(mol);
        CHECK(std::abs(hcorr - hcorr_ans) / hcorr_ans < 5.0e-6);

        const double gibbs_ans = 18.38665762 * cal_to_J * 1000.0;
        double gibbs = Chem::gibbs_energy(mol);
        CHECK(std::abs(gibbs - gibbs_ans) / gibbs_ans < 3.0e-6);
    }

    SECTION("CH2ClCH2Cl")
    {
        Numlib::Vec<double> ans{1.070, 1.745, 30.772}; // C&T (2000) erratum

        std::ifstream from;
        Stdutils::fopen(from, "test_ch2clch2cl.inp");
        Chem::Molecule mol(from);

        Chem::Thermodata td(from);
        auto temp = td.get_temperature();

        for (int i = 0; i < temp.size(); ++i) {
            double res = Chem::qtor(mol, temp(i));
            CHECK(std::abs(res - ans(i)) / ans(i) < 5.0e-2);
        }
    }
}
