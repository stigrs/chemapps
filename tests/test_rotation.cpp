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

TEST_CASE("test_rotation")
{
    using namespace Chem;
    using namespace Numlib;
    using namespace Stdutils;

    std::ifstream from;
    fopen(from, "test_rotation.inp");
    Molecule mol(from);

    Numlib::Vec<double> ans = {127.63201, 24.89071, 24.02767};
    Numlib::Vec<double> res = mol.rot_const();

    for (Index i = 0; i < ans.size(); ++i) {
        CHECK(std::abs(res(i) - ans(i)) < 1.0e-4);
    }
}
