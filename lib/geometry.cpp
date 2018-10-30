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

#include <chem/impl/geometry.h>
#include <stdutils/stdutils.h>

Chem::Impl::Geometry::Geometry(std::istream& from, const std::string& key)
    : atms(), xyz(), zmat(atms, xyz)
{
    using namespace Stdutils;

    auto pos = find_token(from, key);
    if (pos != -1) {
        pos = find_token(from, "geometry", pos);
        if (pos != -1) {
            Impl::read_xyz_format(from, atms, xyz, title);
        }
    }
    if (atms.empty()) {
        pos = find_token(from, key);
        if (pos != -1) {
            from.ignore();
            from.clear();
            zmat.load(from);
        }
    }
}

double Chem::Impl::Geometry::tot_mass() const
{
    double res = 0.0;
    for (const auto& at : atms) {
        res += at.atomic_mass;
    }
    return res;
}

