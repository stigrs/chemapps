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

#ifndef CHEM_ELEMENT_H
#define CHEM_ELEMENT_H

#include <string>

//
// Struct for holding an element in the Periodic Table of Elements.
//
struct Element {
    std::string atomic_symbol;  // atomic or isotopic symbol
    int atomic_number;          // atomic number
    int mass_number;            // mass number
    double atomic_mass;         // atomic mass
    double atomic_weight;       // atomic weight
    double isotope_comp;        // isotopic composition
};

#endif  // CHEM_ELEMENT_H
