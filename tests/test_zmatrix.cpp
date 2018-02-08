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

#include <chem/element.h>
#include <chem/molecule.h>
#include <chem/utils.h>
#include <srs/array.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>
#include <vector>

TEST_CASE("test_zmatrix")
{
    srs::dmatrix xyz_ans = {{0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
                            {1.54000000e+00, 0.00000000e+00, 0.00000000e+00},
                            {2.06671102e+00, 8.86109501e-17, 1.44712664e+00},
                            {2.24685680e+00, 1.44712664e+00, 1.94207310e+00},
                            {-3.72801956e-01, 6.27181400e-17, 1.02426496e+00},
                            {-3.72801956e-01, 8.87039473e-01, -5.12132478e-01},
                            {-3.72801956e-01, -8.87039473e-01, -5.12132478e-01},
                            {1.91280196e+00, -8.87039473e-01, -5.12132478e-01},
                            {1.91280196e+00, 8.87039473e-01, -5.12132478e-01},
                            {3.02776125e+00, -5.12132478e-01, 1.49406052e+00},
                            {1.36067235e+00, -5.12132478e-01, 2.10083125e+00},
                            {1.96127393e+00, 2.14776513e+00, 1.15744062e+00},
                            {1.62151810e+00, 1.62228626e+00, 2.81749906e+00},
                            {3.28860700e+00, 1.62228626e+00, 2.21072833e+00}};

    std::ifstream from;
    srs::fopen(from, "test_zmatrix.inp");

    std::vector<int> moiety = {3, 9, 10};

    Molecule mol(from);
    mol.get_zmat().load(from);
    mol.get_zmat().rotate_moiety(moiety, 90.0);
    srs::dmatrix xyz = mol.get_xyz();

    CHECK(srs::approx_equal(xyz_ans, xyz, 1.0e-8));
}
