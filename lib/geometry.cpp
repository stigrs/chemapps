// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

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

