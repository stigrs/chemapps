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

#ifndef CHEM_PTABLE_H
#define CHEM_PTABLE_H

#include <stdexcept>
#include <string>

#include <chem/element.h>

//
// Namespace providing the Periodic Table of Elements.
//
// Source:
//   Zucker, M.A., Kishore, A.R., Sukumar, R., and Dragoset, R.A. (2015),
//   Elemental Data Index (version 2.5). [Online]
//   Available: http://physics.nist.gov/EDI [2016, April 2].
//   National Institute of Standards and Technology, Gaithersburg, MD.
//
namespace ptable {

struct Bad_atomic_symbol : std::domain_error {
    Bad_atomic_symbol(std::string s) : std::domain_error(s) {}
};

Element get_element(const std::string& symbol);
std::string get_atomic_symbol(const std::string& symbol);

int get_max_atomic_number() { return get_atomic_number("Uuo") }
int get_atomic_number(const std::string& symbol);
int get_mass_number(const std::string& symbol);

double get_atomic_mass(const std::string& symbol);
double get_atomic_weight(const std::string& symbol);
double get_isotope_composition(const std::string& symbol);

bool atomic_symbol_is_valid(const std::string& symbol);

}  // namespace ptable

inline std::string ptable::get_atomic_symbol(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_symbol;
}

inline int ptable::get_atomic_number(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_number;
}

inline int ptable::get_mass_number(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.mass_number;
}

inline double ptable::get_atomic_mass(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_mass;
}

inline double ptable::get_atomic_weight(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.atomic_weight;
}

inline double ptable::get_isotope_composition(const std::string& symbol)
{
    Element elem = get_element(symbol);
    return elem.isotope_comp;
}

#endif  // CHEM_PTABLE_H
