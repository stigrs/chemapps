// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/collision.h>
#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <chem/troe.h>
#include <chem/whirab.h>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <fstream>

TEST_CASE("test_troe")
{
    SECTION("h2o_ar")
    {
        using namespace Numlib::Constants;

        std::ifstream from;
        Stdutils::fopen(from, "test_h2o_ar.inp");

        Chem::Molecule mol(from);

        Numlib::Vec<double> freqs_ans = {3657.05, 1594.59, 3755.79};
        auto freqs = mol.vib().frequencies();
        for (Index i = 0; i < freqs.size(); ++i) {
            CHECK(std::abs(freqs(i) - freqs_ans(i)) < 1.0e-12);
        }

        Chem::Troe troe(from, mol);
        Chem::Collision coll(from);

        double f_anh = troe.f_anharm();
        CHECK(std::abs(f_anh - 2.37) < 1.0e-3);

        double e0 = 118.023 * cal_to_J / icm_to_kJ;
        double a = Chem::Whirab::a_corr(mol, e0);
        CHECK(std::abs(a - 0.989) < 1.0e-2);

        double rho = Chem::Whirab::vibr_density_states(mol, e0);
        rho *= cal_to_J;
        CHECK(std::abs(rho - 16.7) < 5.0e-2);

        double temp = 300.0;
        double kT = R * 1.0e-3 * temp;

        double f_e = troe.f_energy(temp);
        CHECK(std::abs(f_e - 1.01) < 1.0e-3);

        double f_rot = troe.f_rotation(temp);
        CHECK(std::abs(f_rot - 95.0) < 0.8);

        double f_free = troe.f_free_rotor(temp);
        CHECK(std::abs(f_free - 1.0) < 1.0e-12);

        double f_hind = troe.f_hind_rotor(temp);
        CHECK(std::abs(f_hind - 1.0) < 1.0e-12);

        double qvib = Chem::qvib(mol, temp, "V=0");
        CHECK(std::abs(qvib - 1.0) < 5.0e-4);

        double z_lj = coll.lj_coll_freq(temp);
        CHECK((std::abs(z_lj - 1.80e+14) / 1.80e+14) < 1.0e-2);

        rho /= cal_to_J;
        double k0 =
            z_lj * (rho * kT / qvib) * f_anh * f_e * f_rot * f_hind * f_free;
        CHECK((std::abs(k0 - 4.1e+17) / 4.1e+17) < 1.0e-1);

        temp = 2000.0;
        kT = R * 1.0e-3 * temp;

        f_e = troe.f_energy(temp);
        CHECK(std::abs(f_e - 1.06) < 5.0e-3);

        f_rot = troe.f_rotation(temp);
        CHECK(std::abs(f_rot - 14.0) < 0.8);

        f_free = troe.f_free_rotor(temp);
        CHECK(std::abs(f_free - 1.0) < 1.0e-12);

        f_hind = troe.f_hind_rotor(temp);
        CHECK(std::abs(f_hind - 1.0) < 1.0e-12);

        qvib = Chem::qvib(mol, temp, "V=0");
        CHECK(std::abs(qvib - 1.69) < 5.0e-3);

        z_lj = coll.lj_coll_freq(temp);
        CHECK((std::abs(z_lj - 2.92e+14) / 2.92e+14) < 1.0e-3);

        k0 = z_lj * (rho * kT / qvib) * f_anh * f_e * f_rot * f_hind * f_free;
        CHECK((std::abs(k0 - 4.2e+17) / 4.2e+17) < 1.0e-1);
    }
}

