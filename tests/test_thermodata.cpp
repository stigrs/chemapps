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
#include <srs/array.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <string>

TEST_CASE("test_thermodata")
{
    srs::dvector t_ans      = {100.0, 200.0, 300.0, 800.0, 1000.0};
    std::string zeroref_ans = "V=0";

    std::ifstream from;
    srs::fopen(from, "test_thermodata.inp");

    Thermodata td(from);

    CHECK(td.get_vibr_zeroref() == zeroref_ans);
    for (int i = 0; i < t_ans.size(); ++i) {
        CHECK(td.get_temperature()(i) == t_ans(i));
    }
}
