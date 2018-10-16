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

#ifndef CHEM_TRAITS_H
#define CHEM_TRAITS_H

#include <string>

namespace Chem {

// Enumeration of molecular structure types.
enum Mol_type { atom, linear, nonlinear };

// Enumeration of potential types.
enum Pot_type { type1, type2 };

// Struct for holding molecular formula.
struct Mol_formula {
    std::string atom; // atomic or isotopic symbol
    int stoich;       // stoichiometry
};

}  // namespace Chem

#endif  // CHEM_TRAITS_H

