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

#ifndef CHEM_UNITS_H
#define CHEM_UNITS_H

#include <iostream>
#include <string>

namespace Chem {

// Namespace providing methods for handling units.
namespace Units {

enum Type {
    kJ_mol,    // kJ/mol
    kcal_mol,  // kcal/mol
    icm,       // cm**-1
    kelvin,
    hartree,
    hertz,
    eV,
    amu,
    kg,
    au
};

Type lexer(const std::string& unit);
void print(std::ostream& to = std::cout);

}  // namespace Units

} // namespace Chem

#endif // CHEM_UNITS_H

