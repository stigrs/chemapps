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

private:
    // Create Gaussian input file.
    void write_com(const Molecule& mol) const;

    std::string version;   // Gaussian version
    std::string keywords;  // list of Gaussian keywords
    std::string jobname;   // Gaussian job name
    int nprocshared;       // number of processors
    bool nosave;           // flag to specify if chk file should be saved
};

inline Gaussian::Gaussian()
{
    version     = "rung09";
    keywords    = "opt freq hf/sto-3g";
    jobname     = "gauss";
    nprocshared = 1;
    nosave      = true;
}

inline Gaussian::Gaussian(std::istream& from, const std::string& key)
{
    init(from, key);
}

#endif  // CHEM_GAUSSIAN_H
