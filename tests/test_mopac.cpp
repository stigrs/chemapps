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

#include <chem/molecule.h>
#include <chem/mopac.h>
#include <srs/array.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_mopac")
{
    std::ifstream from;
    srs::fopen(from, "test_mopac.inp");
    Molecule mol(from);
    Mopac mop(from);
    mop.run(mol);

    const double heat_ans = 24.49922 * datum::cal_to_J;
    double heat           = mop.get_heat_of_formation();

    CHECK(std::abs(heat - heat_ans) < 2.0e-3);

    srs::dmatrix xyz_ans = {{0.0000, 0.0000, 0.0000},
                            {1.3987, 0.0000, 0.0000},
                            {2.0981, 1.2114, 0.0000},
                            {1.3987, 2.4228, 0.0003},
                            {-0.0001, 2.4229, 0.0008},
                            {-0.6992, 1.2114, -0.0002},
                            {-0.5441, -0.9424, -0.0000},
                            {-0.5441, 3.3653, 0.0011},
                            {1.9428, -0.9424, 0.0000},
                            {-1.7874, 1.2116, -0.0002},
                            {3.1862, 1.2114, 0.0000},
                            {1.9427, 3.3652, 0.0003}};
    srs::dmatrix xyz(xyz_ans.rows(), xyz_ans.cols());
    mop.get_xyz(xyz);
    CHECK(srs::approx_equal(xyz, xyz_ans, 4.0e-3) == true);
}
