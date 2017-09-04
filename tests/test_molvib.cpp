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
#include <chem/element.h>
#include <chem/molecule.h>
#include <chem/utils.h>
#include <armadillo>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_molvib")
{
    double zpe_ans = 0.052023;

    std::ifstream from;
    chem::fopen(from, "test_molvib.inp");

    Molecule mol(from);
    double zpe = mol.get_vib().zero_point_energy() / datum::au_to_icm;

    CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);
}
