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

#ifndef CHEM_WHITTEN_RABINO_H
#define CHEM_WHITTEN_RABINO_H

#include <chem/molecule.h>

//
// Provides Whitten-Rabinovitch approximations.
//
namespace wr {

// Whitten-Rabinovitch correction.
double a_corr(const Molecule& mol, double e_barrier);

// Vibrational density of states.
double vibr_density_states(const Molecule& mol, double e_barrier);

}  // namespace wr

#endif  // CHEM_WHITTEN_RABINO_H
