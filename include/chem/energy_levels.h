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

#ifndef CHEM_ENERGY_LEVELS_H
#define CHEM_ENERGY_LEVELS_H

#include <vector>

namespace Chem {

namespace Energy_levels {

    // Calculate harmonic oscillator energy levels up to a maximum energy.
    std::vector<double> harmonic_oscillator(double freq, double emax);

    // Calculate free rotor energy levels up to a maximum energy.
    std::vector<double> free_rotor(double rotc, double emax);

    // Calculate hindered rotor energy levels up to a maximum energy.
    //
    // For energies less than 1.5 times the barrier height the hindered rotor
    // energy levels are computed using the method of Barker and Shovlin
    // (Chem. Phys. Lett., 2004, vol. 383, pp. 203-207), while free rotor
    // energy levels are used for higher energies and if the barrier height is
    // less than 1.0 cm^-1.
    //
    // Note: Only implemented for one-dimensional hindered rotors.
    //
    std::vector<double>
    hindered_rotor(double sigma, double rotc, double barrier, double emax);

} // namespace Energy_levels

} // namespace Chem

#endif // CHEM_ENERGY_LEVELS_H

