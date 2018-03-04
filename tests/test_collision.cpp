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

#include <chem/collision.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <fstream>
#include <iostream>


TEST_CASE("test_collision")
{
    SECTION("c3h7br_ne_84")
    {
        std::ifstream from;
        srs::fopen(from, "test_collision_c3h7br_ne_84.inp");

        double red_mass    = 20.0 * 121.97 / (20.0 + 121.97);
        double eps_complex = std::sqrt(32.0 * 406.12);
        double sig_complex = 0.5 * (2.82 + 5.78);

        Collision coll(from);

        CHECK(coll.reduced_mass() == red_mass);
        CHECK(coll.average_mass() == 0.0);
        CHECK(coll.epsilon_complex() == eps_complex);
        CHECK(coll.sigma_complex() == sig_complex);
        CHECK(coll.epsilon_local() == 0.0);
        CHECK(coll.sigma_local() == 0.0);
    }

    SECTION("azulene_ne")
    {
        std::ifstream from;
        srs::fopen(from, "test_collision_azulene_ne_90a.inp");

        Collision coll(from);
        coll.biased_random_walk();
    }
}