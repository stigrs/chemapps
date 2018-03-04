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
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <fstream>


TEST_CASE("test_collision")
{
    SECTION("C3H7Br-Ne")
    {
        std::ifstream from;
        srs::fopen(from, "test_collision_c3h7br_ne.inp");

        Collision coll(from);
        coll.biased_random_walk();
        CHECK(srs::approx_equal(coll.sigma_complex(), 4.30, 1.0e-12));
        CHECK(srs::approx_equal(coll.epsilon_complex(), 113.999, 1.0e-3));
        CHECK(srs::approx_equal(coll.reduced_mass(), 17.183, 1.0e-3));
    }
}