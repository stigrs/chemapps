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

#include <chem/collision.h>
#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <chem/troe.h>
#include <chem/whitten_rabino.h>
#include <srs/array.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <cmath>
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
        Collision coll(from);

        double f_anh = troe.f_anharm();
        CHECK(srs::approx_equal(f_anh, 2.37, 1.0e-3));

        double e0 = 118.023 * datum::cal_to_J / datum::icm_to_kJ;
        double a  = wr::a_corr(mol, e0);
        CHECK(srs::approx_equal(a, 0.989, 1.0e-2));

        double rho = wr::vibr_density_states(mol, e0);
        rho *= datum::cal_to_J;
        CHECK(srs::approx_equal(rho, 16.7, 5.0e-2));

        double temp = 300.0;

        double kT = datum::R * 1.0e-3 * temp;

        double f_e = troe.f_energy(temp);
        CHECK(srs::approx_equal(f_e, 1.01, 1.0e-3));

        double f_rot = troe.f_rotation(temp);
        CHECK(srs::approx_equal(f_rot, 95.0, 0.8));

        double f_free = troe.f_free_rotor(temp);
        CHECK(srs::approx_equal(f_free, 1.0, 1.0e-12));

        double f_hind = troe.f_hind_rotor(temp);
        CHECK(srs::approx_equal(f_hind, 1.0, 1.0e-12));

        double qvib = chem::qvib(mol, temp, "V=0");
        CHECK(srs::approx_equal(qvib, 1.0, 5.0e-4));

        double z_lj = coll.lj_coll_freq(temp);
        CHECK(srs::approx_equal(z_lj, 1.80e+14, 1.0e-2, "reldiff"));

        rho /= datum::cal_to_J;
        double k0
            = z_lj * (rho * kT / qvib) * f_anh * f_e * f_rot * f_hind * f_free;
        CHECK(srs::approx_equal(k0, 4.1e+17, 1.0e-1, "reldiff"));

        temp = 2000.0;

        kT = datum::R * 1.0e-3 * temp;

        f_e = troe.f_energy(temp);
        CHECK(srs::approx_equal(f_e, 1.06, 5.0e-3));

        f_rot = troe.f_rotation(temp);
        CHECK(srs::approx_equal(f_rot, 14.0, 0.8));

        f_free = troe.f_free_rotor(temp);
        CHECK(srs::approx_equal(f_free, 1.0, 1.0e-12));

        f_hind = troe.f_hind_rotor(temp);
        CHECK(srs::approx_equal(f_hind, 1.0, 1.0e-12));

        qvib = chem::qvib(mol, temp, "V=0");
        CHECK(srs::approx_equal(qvib, 1.69, 5.0e-3));

        z_lj = coll.lj_coll_freq(temp);
        CHECK(srs::approx_equal(z_lj, 2.92e+14, 1.0e-3, "reldiff"));

        k0 = z_lj * (rho * kT / qvib) * f_anh * f_e * f_rot * f_hind * f_free;
        CHECK(srs::approx_equal(k0, 4.2e+17, 1.0e-1, "reldiff"));
    }
}
