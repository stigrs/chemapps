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
#include <srs/utils.h>
#include <map>

Thermodata::Thermodata(std::istream& from, const std::string& key)
{
    // Read input data:

    srs::dvector p_def = {datum::std_atm};
    srs::dvector t_def = {298.15};

    std::map<std::string, srs::Input> input_data;
    input_data["pressure"]    = srs::Input(pressure, p_def);
    input_data["temperature"] = srs::Input(temperature, t_def);
    input_data["incl_sigma"]  = srs::Input(incl_sigma, 1);
    input_data["zeroref"]     = srs::Input(zeroref, "BOT");

    bool found = srs::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }

    // Check if initialized:

    for (auto& it : input_data) {
        if (!it.second.is_init()) {
            throw Thermodata_error(it.first + " not initialized");
        }
    }
    // TODO (stigrs@gmail.com) Implement validation
}
