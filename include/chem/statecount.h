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

#ifndef CHEM_STATECOUNT_H
#define CHEM_STATECOUNT_H

#include <chem/molecule.h>
#include <numlib/matrix.h>

namespace Chem {

namespace Statecount {

    // Count density or sum of states for a molecule.
    //
    // Args:
    //   mol: molecule object
    //   ngrains: the number of energy grains
    //   egrain: energy grain size (cm^-1)
    //   sum: flag to specify if sum of states should be computed
    //
    // Return:
    //   array with rovibrational density or sum of states
    //
    Numlib::Vec<double> count(const Molecule& mol,
                              int ngrains,
                              double egrain = 1.0,
                              bool sum = false);

    // Modified Beyer-Swinehart algorithm for the rovibrational density or sum
    // of states of a system of n harmonic oscillators.
    //
    // Algorithm:
    //   Tables on page 157 and 158 in Gilbert and Smith (1990).
    //
    // Args:
    //   vibr: collection of n harmonic oscillators
    //   ngrains: the number of energy grains
    //   egrain: energy grain size (cm^-1)
    //   sum: flag to specify if sum of states should be computed
    //   rot: array with density or sum of rotational states
    //
    // Returns:
    //   rovibrational density or sum of states
    //
    Numlib::Vec<double>
    bswine(const Numlib::Vec<double>& vibr,
           int ngrains,
           double egrain = 1.0,
           bool sum = false,
           const Numlib::Vec<double>& rot = Numlib::Vec<double>{});

    Numlib::Vec<double> steinrab(const Numlib::Vec<double>& vibr,
                                 double sigma,
                                 double rotc,
                                 double barrier,
                                 int ngrains,
                                 double egrain = 1.0,
                                 bool sum = false);

    // Calculate the density or sum of states for one independent free rotor.
    //
    // Algorithm:
    //   Eq. 4.19 in Forst (2003)
    //
    // Args:
    //   sigma: symmetry number for the free rotor
    //   rotc: rotational constant (cm^-1)
    //   ngrains: the number of energy grains
    //   egrain: energy grain size (cm^-1)
    //   sum: flag to specify if sum of states should be computed
    //
    // Returns:
    //   density or sum of states for the free rotor
    //
    Numlib::Vec<double> free_rotor(double sigma,
                                   double rotc,
                                   int ngrains,
                                   double egrain = 1.0,
                                   bool sum = false);

    // Calculate the density or sum of states for a classical 1D hindered rotor.
    //
    // Algorithm:
    //   Eqs. 4.52 and 4.53 in Forst (2003)
    //
    // Args:
    //   sigma: symmetry number for the free rotor
    //   rotc: rotational constant (cm^-1)
    //   barrier: torsional barrier (cm^-1)
    //   ngrains: the number of energy grains
    //   egrain: energy grain size (cm^-1)
    //   sum: flag to specify if sum of states should be computed
    //
    // Returns:
    //   density or sum of states for the hindered rotor
    //
    Numlib::Vec<double> hindered_rotor(double sigma,
                                       double rotc,
                                       double barrier,
                                       int ngrains,
                                       double egrain = 1.0,
                                       bool sum = false);

} // namespace Statecount

} // namespace Chem

#endif // CHEM_STATECOUNT_H

