////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009-2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/units.h>


Units::Type Units::lexer(const std::string& unit)
{
    enum Type ans;
    if ((unit == "kJ/mol") || (unit == "kJ mol**-1")) {
        ans = kJ_mol;
    }
    else if ((unit == "kcal/mol") || (unit == "kcal mol**-1")) {
        ans = kcal_mol;
    }
    else if ((unit == "cm**-1") || (unit == "cm^-1") || (unit == "cm-1")) {
        ans = icm;
    }
    else if ((unit == "kelvin") || (unit == "K")) {
        ans = kelvin;
    }
    else if ((unit == "hartree") || (unit == "Eh")) {
        ans = hartree;
    }
    else if ((unit == "hertz") || (unit == "Hertz") || (unit == "s**-1")
             || (unit == "s^-1") || (unit == "s-1")) {
        ans = hertz;
    }
    else if (unit == "eV") {
        ans = eV;
    }
    else if (unit == "amu") {
        ans = amu;
    }
    else if (unit == "kg") {
        ans = kg;
    }
    else if ((unit == "au") || (unit == "a.u.")) {
        ans = au;
    }
    else {
        throw Unit_error("unknown unit: " + unit);
    }
    return ans;
}

void Units::print(std::ostream& to)
{
    to << "Supported units:\n"
       << " kJ/mol, kJ mol**-1\n"
       << " kcal/mol, kcal mol**-1\n"
       << " cm**-1, cm^-1, cm-1\n"
       << " kelvin, K\n"
       << " hartree, Eh\n"
       << " hertz, Hertz, s**-1, s^-1, s-1\n"
       << " eV\n"
       << " amu\n"
       << " kg\n"
       << " au, a.u.\n";
}
