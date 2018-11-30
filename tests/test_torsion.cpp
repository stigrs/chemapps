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
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_torsion")
{
    using namespace Chem;
    using namespace Stdutils;

    SECTION("CH2ClCH2Cl")
    {
        std::ifstream from;
        fopen(from, "test_ch2clch2cl.inp");

        Molecule mol(from);
        double rmi = mol.tor().red_moment();

        const double rmi_ans = 58.76991427; // Chuang and Truhlar (2000)
        CHECK(std::abs(rmi - rmi_ans) < 9.0e-2);
    }

    SECTION("CH3OH")
    {
        std::ifstream from;
        fopen(from, "test_ch3oh.inp");

        Molecule mol(from);
        double rmi = mol.tor().red_moment();

        const double rmi_ans = 2.19; // Chuang and Truhlar (2000)
        CHECK(std::abs(rmi - rmi_ans) < 2.0e-2);
    }
}
