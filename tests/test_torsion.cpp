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

#include <chem/molecule.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_torsion")
{
    SECTION("CH2ClCH2Cl")
    {
        std::ifstream from;
        srs::fopen(from, "test_torsion_ch2clch2cl.inp");

        Molecule mol(from);
        double rmi = mol.get_tor().red_moment_of_inertia();

        const double rmi_ans = 58.76991427;  // Chuang and Truhlar (2000)
        CHECK(std::abs(rmi - rmi_ans) < 9.0e-2);
    }

    SECTION("CH3OH")
    {
        std::ifstream from;
        srs::fopen(from, "test_torsion_ch3oh.inp");

        Molecule mol(from);
        double rmi = mol.get_tor().red_moment_of_inertia();

        const double rmi_ans = 2.19;  // Chuang and Truhlar (2000)
        CHECK(std::abs(rmi - rmi_ans) < 2.0e-2);
    }
}
