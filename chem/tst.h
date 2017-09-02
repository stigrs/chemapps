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

#ifndef CHEM_TST_H
#define CHEM_TST_H

#include <chem/molecule.h>
#include <chem/tunnel.h>
#include <iostream>
#include <memory>
#include <string>

//
// Class providing Transition State Theory (TST).
//
// Note: Currently, only conventional TST is implemented. It is recommended
// to use Polyrate if variational TST is needed.
//
class Tst {
public:
    Tst(std::istream& from, const std::string& key = "TST");

    ~Tst() {}

private:
    enum Method_t { conventional };
    enum Reaction_t { unimolecular, bimolecular };

    Method_t method     = conventional;  // TST method
    Reaction_t reaction = bimolecular;   // reaction type

    std::unique_ptr<Molecule> ra;  // reactant A
    std::unique_ptr<Molecule> rb;  // reactant B
    std::unique_ptr<Molecule> ts;  // transition state

    double en_barrier;  // reaction barrier (kJ/mol)
    int rxn_sigma;      // reaction symmetry number
};

#endif  // CHEM_TST_H
