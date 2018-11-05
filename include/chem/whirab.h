////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#ifndef CHEM_WHIRAB_H
#define CHEM_WHIRAB_H

#include <chem/molecule.h>

namespace Chem {

// Provides Whitten-Rabinovitch approximations.
//
namespace Whirab {

// Whitten-Rabinovitch correction.
double a_corr(const Molecule& mol, double e_barrier);

// Vibrational density of states.
double vibr_density_states(const Molecule& mol, double e_barrier);

} // namespace Whirab

} // namespace Chem

#endif // CHEM_WHIRAB_H

