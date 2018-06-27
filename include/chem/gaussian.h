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

#ifndef CHEM_GAUSSIAN_H
#define CHEM_GAUSSIAN_H

#include <chem/gauss_data.h>
#include <chem/molecule.h>
#include <iostream>
#include <string>

//
// Wrapper class for running Gaussian calculations.
//
class Gaussian {
public:
    Gaussian();

    Gaussian(std::istream& from, const std::string& key = "Gaussian");

    // Initialize Gaussian calculation.
    void init(std::istream& from, const std::string& key = "Gaussian");

    // Run Gaussian calculation.
    void run(Molecule& mol) const;
};

#endif  // CHEM_GAUSSIAN_H
