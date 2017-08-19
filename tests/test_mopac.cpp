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

#include <chem/datum.h>
#include <chem/molecule.h>
#include <chem/mopac.h>
#include <chem/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_mopac")
{
    std::ifstream from;
    chem::fopen(from, "test_mopac.inp");
    Molecule mol(from);
    Mopac mop(from);
    mop.run(mol);

    const double heat_ans = 112.00281 * datum::cal_to_J;
    double heat           = mop.get_heat_of_formation();

    CHECK(std::abs(heat - heat_ans) < 1.0e-12);
}
