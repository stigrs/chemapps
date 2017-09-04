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

#include <chem/thermodata.h>
#include <chem/tst.h>
#include <chem/utils.h>
#include <armadillo>
#include <catch/catch.hpp>
#include <fstream>

TEST_CASE("test_tst")
{
    SECTION("ch4oh")
    {
        std::ifstream from;
        chem::fopen(from, "test_tst_ch4oh.inp");

        Tst tst(from);
        Thermodata td(from);

        arma::vec temp = td.get_temperature();

        tst.rate();
        for (arma::uword i = 0; i < temp.size(); ++i) {
            std::cout << tst.rate_coeff(temp(i)) << '\n';
        }
    }
}