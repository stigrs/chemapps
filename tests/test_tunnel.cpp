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
#include <chem/tunnel.h>
#include <chem/utils.h>
#include <armadillo>
#include <catch/catch.hpp>
#include <fstream>

TEST_CASE("test_tunnel")
{
    SECTION("eckart")
    {
        std::ifstream from;
        chem::fopen(from, "test_tunnel.inp");

        Tunnel tunnel(from);
        Thermodata td(from);

        // Results from Brown (1981):
        arma::vec kappa_ans = {6.3946, 2.4444, 1.3019, 1.1761, 1.1189};
        arma::vec temp      = td.get_temperature();

        for (arma::uword i = 0; i < temp.size(); ++i) {
            double kappa = tunnel.eckart(temp(i));
            CHECK(std::abs(kappa - kappa_ans(i)) < 5.0e-3);
        }
    }
}