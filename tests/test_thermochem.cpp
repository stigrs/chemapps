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

TEST_CASE("test_thermochem")
{
    std::ifstream from;
    chem::fopen(from, "test_thermochem_ch3oh.inp");
    Molecule mol(from);

    SECTION("CH3OH")
    {
        const double qtr_ans = 0.712383e+7;
        double qtr           = chem::qtrans(mol, 298.15, datum::std_atm);
        CHECK(std::abs(qtr - qtr_ans) / qtr_ans < 5.0e-6);

        const double str_ans = 36.324 * datum::cal_to_J;
        double str           = chem::entropy_trans(mol);
        CHECK(std::abs(str - str_ans) / str_ans < 1.0e-8);
    }
}