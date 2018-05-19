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

#ifndef CHEM_THERMODATA_H
#define CHEM_THERMODATA_H

#include <srs/array.h>
#include <srs/datum.h>
#include <iostream>
#include <stdexcept>
#include <string>

// Error reporting:

struct Thermodata_error : std::runtime_error {
    Thermodata_error(std::string s) : std::runtime_error(s) {}
};
//
// Class for holding thermochemical data.
//
class Thermodata {
public:
    Thermodata();
    Thermodata(std::istream& from, const std::string& key = "ThermoData");

    srs::dvector get_pressure() const { return pressure; }
    srs::dvector get_temperature() const { return temperature; }
    bool incl_rot_symmetry() const { return incl_sigma; }
    std::string get_vibr_zeroref() const { return zeroref; }

    void set_pressure(const srs::dvector& p) { pressure = p; }
    void set_temperature(const srs::dvector& t) { temperature = t; }
    void set_incl_rot_symmetry(bool flag) { incl_sigma = flag; }
    void set_vibr_zeroref(const std::string& ref) { zeroref = ref; }

private:
    srs::dvector pressure;
    srs::dvector temperature;
    int incl_sigma;       // include rotational symmetry number?
    std::string zeroref;  // zero reference point for vibrations (BOT or V=0)
};

inline Thermodata::Thermodata()
{
    pressure    = {datum::std_atm};
    temperature = {298.15};
    incl_sigma  = 1;
    zeroref     = "BOT";
}

#endif  // CHEM_THERMODATA_H
