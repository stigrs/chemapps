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
#include <armadillo>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_molrot")
{
    arma::vec3 rotc_ans = {127.63201, 24.89071, 24.02767};

    std::ifstream from;
    chem::fopen(from, "test_molrot.inp");

    Molecule mol(from);
    arma::vec3 rotc = mol.get_rot().constants();

    for (arma::uword i = 0; i < rotc.size(); ++i) {
        CHECK(std::abs(rotc(i) - rotc_ans(i)) < 1.0e-4);
    }
}
