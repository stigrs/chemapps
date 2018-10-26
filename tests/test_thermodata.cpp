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
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <string>
#include <cmath>

TEST_CASE("test_thermodata")
{
    Numlib::Vec<double> t_ans = {100.0, 200.0, 300.0, 800.0, 1000.0};
    std::string zeroref_ans = "V=0";

    std::ifstream from;
    Stdutils::fopen(from, "test_thermodata.inp");

    Chem::Thermodata td(from);

    CHECK(td.get_vibr_zeroref() == zeroref_ans);
    for (Index i = 0; i < t_ans.size(); ++i) {
        CHECK(std::abs(td.get_temperature()(i) - t_ans(i)) < 1.0e-12);
    }
}

