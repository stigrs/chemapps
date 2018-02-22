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

#include <chem/element.h>
#include <chem/molecule_io.h>
#include <srs/array.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <fstream>
#include <string>
#include <vector>

TEST_CASE("test_molecule_io")
{
    srs::dmatrix xyz_ans(8, 3);
    xyz_ans(0, 0) = 0.0000;
    xyz_ans(0, 1) = 0.0000;
    xyz_ans(0, 2) = 0.7637;
    xyz_ans(1, 0) = 0.0000;
    xyz_ans(1, 1) = 0.0000;
    xyz_ans(1, 2) = -0.7637;
    xyz_ans(2, 0) = 0.0000;
    xyz_ans(2, 1) = 1.0121;
    xyz_ans(2, 2) = 1.1564;
    xyz_ans(3, 0) = -0.8765;
    xyz_ans(3, 1) = -0.5060;
    xyz_ans(3, 2) = 1.1564;
    xyz_ans(4, 0) = 0.8765;
    xyz_ans(4, 1) = -0.5060;
    xyz_ans(4, 2) = 1.1564;
    xyz_ans(5, 0) = 0.0000;
    xyz_ans(5, 1) = -1.0121;
    xyz_ans(5, 2) = -1.1564;
    xyz_ans(6, 0) = -0.8765;
    xyz_ans(6, 1) = 0.5060;
    xyz_ans(6, 2) = -1.1564;
    xyz_ans(7, 0) = 0.8765;
    xyz_ans(7, 1) = 0.5060;
    xyz_ans(7, 2) = -1.1564;

    std::string title;
    srs::dmatrix xyz;
    std::vector<Element> atoms;

    std::ifstream from;
    srs::fopen(from, "test_molecule_io.inp");
    chem::read_xyz_format(from, atoms, xyz, title);

    CHECK(srs::approx_equal(xyz, xyz_ans, 1.0e-12));
}
