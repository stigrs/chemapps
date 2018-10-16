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
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_molecule")
{
    using namespace Chem;
    using namespace Numlib;
    using namespace Stdutils;

    std::ifstream from;
    fopen(from, "test_molecule.inp");
    Molecule mol(from, "Molecule", false);

    SECTION("num_atoms") { CHECK(mol.num_atoms() == 8); }

    SECTION("net_charge") { CHECK(mol.net_charge() == 0); }

    SECTION("spin_mult") { CHECK(mol.spin_mult() == 2); }

    SECTION("elec_energy")
    {
        CHECK(std::abs(mol.elec_energy() + 0.5) < 1.0e-12);
    }

    SECTION("spin-orbit")
    {
        Vec<int> degen_ans = {2, 2};
        Vec<double> energy_ans = {0.0, 140.0};

        CHECK(same_extents(mol.spin_orbit_degen(), mol.spin_orbit_energy()));

        CHECK(mol.spin_orbit_degen() == degen_ans);

        for (Index i = 0; i < degen_ans.size(); ++i) {
            CHECK(std::abs(mol.spin_orbit_energy()(i) - energy_ans(i)) <
                  1.0e-12);
        }
    }

    SECTION("geometry")
    {
        Mat<double> ans = {
            {0.0000, 0.0000, 0.7637},   {0.0000, 0.0000, -0.7637},
            {0.0000, 1.0121, 1.1564},   {-0.8765, -0.5060, 1.1564},
            {0.8765, -0.5060, 1.1564},  {0.0000, -1.0121, -1.1564},
            {-0.8765, 0.5060, -1.1564}, {0.8765, 0.5060, -1.1564}};

        auto res = mol.cart_coord();
        CHECK(same_extents(res, ans));

        for (Index i = 0; i < ans.rows(); ++i) {
            for (Index j = 0; j < ans.cols(); ++j) {
                CHECK(std::abs(res(i, j) - ans(i, j)) < 1.0e-12);
            }
        }
    }
}
