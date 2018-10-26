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
#include <stdutils/stdutils.h>

Chem::Thermodata::Thermodata(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    // Read input data:

    Numlib::Vec<double> p_def = {Numlib::Constants::std_atm};
    Numlib::Vec<double> t_def = {298.15};

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "pressure", pressure, p_def);
        get_token_value(from, pos, "temperature", temperature, t_def);
        get_token_value(from, pos, "incl_sigma", incl_sigma, 1);
        get_token_value(from, pos, "zeroref", zeroref, std::string("BOT"));
    }

    // Validate input:

    for (const auto& p : pressure) {
        Assert::dynamic(p > 0.0, "bad pressure");
    }
    for (const auto& t : temperature) {
        Assert::dynamic(t > 0.0, "bad temperature");
    }
    Assert::dynamic(incl_sigma == 0 || incl_sigma == 1, "bad incl_sigma");
    Assert::dynamic(zeroref == "BOT" || zeroref == "V=0", "bad zeroref");
}

