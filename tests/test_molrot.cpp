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
#include <srs/utils.h>
#include <srs/array.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_molrot")
{
    srs::dvector rotc_ans = {127.63201, 24.89071, 24.02767};

    std::ifstream from;
    srs::fopen(from, "test_molrot.inp");

    Molecule mol(from);
    srs::dvector rotc = mol.get_rot().constants();

    for (int i = 0; i < rotc.size(); ++i) {
        CHECK(std::abs(rotc(i) - rotc_ans(i)) < 1.0e-4);
    }
}
