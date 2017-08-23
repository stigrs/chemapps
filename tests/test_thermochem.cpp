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
#include <chem/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

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
    }
}