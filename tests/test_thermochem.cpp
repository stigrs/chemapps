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
#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <chem/thermodata.h>
#include <chem/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_thermochem")
{
    SECTION("CH3OH")
    {
        std::ifstream from;
        chem::fopen(from, "test_thermochem_ch3oh.inp");
        Molecule mol(from);

        const double qtr_ans = 0.712383e+7;
        double qtr           = chem::qtrans(mol, 298.15, datum::std_atm);
        CHECK(std::abs(qtr - qtr_ans) / qtr_ans < 5.0e-6);

        const double htr_ans = 0.889 * datum::cal_to_J * 1000.0;
        double htr           = chem::thermal_energy_trans();
        CHECK(std::abs(htr - htr_ans) / htr_ans < 5.0e-4);

        const double str_ans = 36.324 * datum::cal_to_J;
        double str           = chem::entropy_trans(mol);
        CHECK(std::abs(str - str_ans) / str_ans < 1.0e-8);

        const double cv_tr_ans = 2.981 * datum::cal_to_J;
        double cv_tr           = chem::const_vol_heat_trans();
        CHECK(std::abs(cv_tr - cv_tr_ans) / cv_tr_ans < 1.0e-4);

        const double qe_ans = 0.100000e+1;
        double qe           = chem::qelec(mol);
        CHECK(std::abs(qe - qe_ans) / qe_ans < 1.0e-8);

        const double he_ans = 0.0;
        double he           = chem::thermal_energy_elec();
        CHECK(std::abs(he - he_ans) < 1.0e-9);

        const double cv_e_ans = 0.0;
        double cv_e           = chem::const_vol_heat_elec();
        CHECK(std::abs(cv_e - cv_e_ans) < 1.0e-9);

        const double qr_ans = 0.314134e+4;
        double qr           = chem::qrot(mol);
        CHECK(std::abs(qr - qr_ans) / qr_ans < 1.0e-5);

        const double sr_ans = 18.983 * datum::cal_to_J;
        double sr           = chem::entropy_rot(mol);
        CHECK(std::abs(sr - sr_ans) / sr_ans < 3.0e-5);

        const double er_ans = 0.889 * datum::cal_to_J * 1000.0;
        double er           = chem::thermal_energy_rot(mol);
        CHECK(std::abs(er - er_ans) / er_ans < 5.0e-4);

        const double cv_r_ans = 2.981 * datum::cal_to_J;
        double cv_r           = chem::const_vol_heat_rot(mol);
        CHECK(std::abs(cv_r - cv_r_ans) / cv_r_ans < 1.0e-4);

        const double qv_ans = 0.148820e-23;
        double qv           = chem::qvib(mol);
        CHECK(std::abs(qv - qv_ans) / qv_ans < 1.0e-4);

        const double sv_ans = 1.441 * datum::cal_to_J;
        double sv           = chem::entropy_vib(mol);
        CHECK(std::abs(sv - sv_ans) / sv_ans < 2.0e-4);

        const double ev_ans = 32.936 * datum::cal_to_J * 1000.0;
        double ev           = chem::thermal_energy_vib(mol);
        CHECK(std::abs(ev - ev_ans) / ev_ans < 5.0e-5);

        const double cv_v_ans = 2.696 * datum::cal_to_J;
        double cv_v           = chem::const_vol_heat_vib(mol);
        CHECK(std::abs(cv_v - cv_v_ans) / cv_v_ans < 2.0e-4);

        const double qtot_ans = 0.333036e-13;
        double qtot           = chem::qtot(mol);
        CHECK(std::abs(qtot - qtot_ans) / qtot_ans < 1.0e-4);

        const double stot_ans = 56.748 * datum::cal_to_J;
        double stot           = chem::entropy(mol);
        CHECK(std::abs(stot - stot_ans) / stot_ans < 1.0e-5);

        const double etot_ans = 34.714 * datum::cal_to_J * 1000.0;
        double etot           = chem::thermal_energy(mol);
        CHECK(std::abs(etot - etot_ans) / etot_ans < 2.0e-5);

        const double cv_tot_ans = 8.657 * datum::cal_to_J;
        double cv_tot           = chem::const_vol_heat_capacity(mol);
        CHECK(std::abs(cv_tot - cv_tot_ans) / cv_tot_ans < 3.0e-5);

        const double hcorr_ans = 35.3061979 * datum::cal_to_J * 1000.0;
        double hcorr           = chem::enthalpy(mol);
        CHECK(std::abs(hcorr - hcorr_ans) / hcorr_ans < 5.0e-6);

        const double gibbs_ans = 18.38665762 * datum::cal_to_J * 1000.0;
        double gibbs           = chem::gibbs_energy(mol);
        CHECK(std::abs(gibbs - gibbs_ans) / gibbs_ans < 3.0e-6);
    }

    SECTION("CH2ClCH2Cl")
    {
        arma::vec qtor_ans{1.070, 1.745, 30.772};  // C&T (2000) erratum

        std::ifstream from;
        chem::fopen(from, "test_thermochem_ch2clch2cl.inp");
        Molecule mol(from);

        Thermodata td(from);
        arma::vec temp = td.get_temperature();

        for (arma::uword i = 0; i < temp.size(); ++i) {
            double qtor = chem::qtor(mol, temp(i));
            CHECK(std::abs(qtor - qtor_ans(i)) / qtor_ans(i) < 5.0e-2);
        }
    }
}
