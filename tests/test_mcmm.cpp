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

#include <chem/mcmm.h>
#include <chem/molecule.h>
#include <chem/mopac.h>
#include <srs/datum.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>


TEST_CASE("test_mcmm")
{
    std::ifstream from;
    srs::fopen(from, "test_mcmm.inp");

    Molecule mol(from);
    mol.get_zmat().load(from);
    Mcmm<Mopac> mc(from, mol);
    double eglobal   = mc.get_global_min_energy();
    double emin_anti = -31.03257 * datum::cal_to_J;
    CHECK(std::abs(eglobal - emin_anti) < 2.0e-4);
}
