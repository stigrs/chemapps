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

#include <chem/thermodata.h>
#include <chem/tunnel.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <fstream>

TEST_CASE("test_tunnel")
{
    SECTION("eckart")
    {
        std::ifstream from;
        Stdutils::fopen(from, "test_tunnel.inp");

        Chem::Tunnel tunnel(from);
        Chem::Thermodata td(from);

        // Results from Brown (1981):
        Numlib::Vec<double> ans = {6.41025, 2.42735, 1.29311, 1.17025, 1.11455};
        auto temp = td.get_temperature();

        for (int i = 0; i < temp.size(); ++i) {
            double res = tunnel.factor(temp(i));
            CHECK(std::abs(res - ans(i)) < 1.0e-3);
        }
    }

    SECTION("none")
    {
        Chem::Tunnel tunnel_none;
        CHECK(tunnel_none.factor() == 1.0);
    }
}

