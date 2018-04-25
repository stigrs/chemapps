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

#include <chem/ptable.h>
#include <catch/catch.hpp>

TEST_CASE("test_ptable")
{
    SECTION("Carbon")
    {
        double atomic_weight = ptable::get_element("C").atomic_weight;
        CHECK(atomic_weight == (12.0096 + 12.0116) * 0.5);
    }

    SECTION("109Ag")
    {
        double atomic_mass = ptable::get_element("109Ag").atomic_mass;
        CHECK(atomic_mass == 108.904755);
    }

    SECTION("atomic_number") { CHECK(ptable::get_atomic_symbol(15) == "P"); }
}
