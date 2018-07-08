////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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
#include <chem/troe.h>
#include <chem/whitten_rabino.h>
#include <srs/array.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <fstream>


TEST_CASE("test_troe")
{
    SECTION("h2o_ar")
    {
        std::ifstream from;
        srs::fopen(from, "test_troe_h2o_ar.inp");

        Molecule mol(from);

        srs::dvector freqs_ans = {3657.05, 1594.59, 3755.79};
        auto freqs             = mol.get_vib().get_freqs();
        CHECK(srs::approx_equal(freqs, freqs_ans, 1.0e-12));

        Troe troe(from, mol);

        double f_anh = troe.f_anharm();
        CHECK(srs::approx_equal(f_anh, 2.37, 1.0e-3));

        double e0 = 118.023 * datum::cal_to_J / datum::icm_to_kJ;
        double a  = wr::a_corr(mol, e0);
        CHECK(srs::approx_equal(a, 0.989, 1.0e-2));
    }
}
