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

#include <chem/impl/elec_state.h>
#include <stdutils/stdutils.h>
#include <cassert>

Chem::Impl::Elec_state::Elec_state(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    so_degen = {1};
    so_energy = {0.0};

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "charge", charge, 0);
        get_token_value(from, pos, "spin_mult", spin, 1);
        get_token_value(from, pos, "elec_energy", energy, 0.0);
        get_token_value(from, pos, "so_degen", so_degen, so_degen);
        get_token_value(from, pos, "so_energy", so_energy, so_energy);
    }
    assert(same_extents(so_degen, so_energy));
}

