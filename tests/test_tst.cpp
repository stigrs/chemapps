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

#include <chem/tst.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>


TEST_CASE("test_tst")
{
    SECTION("ch4cl")
    {
        std::ifstream from;
        Stdutils::fopen(from, "test_tst_ch4cl.inp");

        // These values are from Polyrate 2017:
        Numlib::Vec<double> ktst_ans = {1.9919E-15, 4.5171E-14, 1.9234E-12,
                                        6.4847E-11};
        Numlib::Vec<double> ktstw_ans = {5.8627E-15, 8.4186E-14, 2.3387E-12,
                                         6.7087E-11};

        Chem::Tst tst(from);
        Chem::Thermodata td(from);

        auto temp = td.get_temperature();

        for (Index i = 0; i < temp.size(); ++i) {
            double ktst  = tst.rate_coeff(temp(i));
            double ktstw = ktst * tst.tunneling(temp(i));
            CHECK(std::abs(ktst - ktst_ans(i)) / ktst_ans(i) < 1.0e-4);
            CHECK(std::abs(ktstw - ktstw_ans(i)) / ktstw_ans(i) < 1.0e-4);
        }
    }

    SECTION("ch4oh")
    {
        std::ifstream from;
        Stdutils::fopen(from, "test_tst_ch4oh.inp");

        // These values are from Polyrate 2017:
        Numlib::Vec<double> ktst_ans = {3.9039E-18, 5.1396E-16, 1.1748E-13,
                                        1.2988E-11, 8.5959E-11};
        Numlib::Vec<double> ktstw_ans = {3.2091E-17, 2.1633E-15, 2.1172E-13,
                                         1.4655E-11, 9.0269E-11};

        Chem::Tst tst(from);
        Chem::Thermodata td(from);

        auto temp = td.get_temperature();

        for (Index i = 0; i < temp.size(); ++i) {
            double ktst  = tst.rate_coeff(temp(i));
            double ktstw = ktst * tst.tunneling(temp(i));
            CHECK(std::abs(ktst - ktst_ans(i)) / ktst_ans(i) < 5.0e-4);
            CHECK(std::abs(ktstw - ktstw_ans(i)) / ktstw_ans(i) < 5.0e-4);
        }
    }
}

