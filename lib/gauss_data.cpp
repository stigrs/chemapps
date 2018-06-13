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

#include <chem/gauss_data.h>
#include <srs/utils.h>
#include <sstream>

Gauss_version Gauss_data::get_version() const
{
    const char pattern_start[] = "Cite this work as:";
    const char pattern_ver[]   = "Gaussian";

    std::streampos orig_pos = from.tellg();

    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file
    from.clear();

    Gauss_version version = unknown;

    std::string line;
    std::string word;
    std::string ver;

    while (std::getline(from, line)) {
        if (line.find(pattern_start, 0) != std::string::npos) {
            std::getline(from, line);
            std::istringstream iss(line);
            iss >> word >> ver;
            if (word == pattern_ver) {
                if (ver == "94,") {
                    version = g94;
                }
                else if (ver == "98,") {
                    version = g98;
                }
                else if (ver == "03,") {
                    version = g03;
                }
                else if (ver == "09,") {
                    version = g09;
                }
            }
        }
    }
    if (version == unknown) {
        throw Gauss_error("unknown Gaussian version");
    }
    from.clear();
    from.seekg(orig_pos, std::ios_base::beg);  // move to original position
    from.clear();

    return version;
}

int Gauss_data::get_natoms() const
{
    std::streampos orig_pos = from.tellg();

    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file
    from.clear();

    int natoms = 0;

    if (filetype == out) {
        const char pattern_start[] = "Z-Matrix";
        const char pattern_end[]   = "Distance";

        std::string line;
        std::string token;

        while (std::getline(from, line)) {
            std::istringstream iss1(line);
            iss1 >> token;
            if (token == pattern_start) {
                for (int i = 0; i < 4; ++i) {  // ignore four lines
                    from.ignore(256, '\n');
                }
                while (std::getline(from, line)) {
                    std::istringstream iss2(line);
                    if (line[1] == '-') {
                        break;
                    }
                    else {
                        iss2 >> natoms;
                    }
                }
            }
            else if (token == pattern_end) {
                break;
            }
        }
    }
    else {  // filetype == fchk
        const char pattern[] = "Number of atoms";

        std::string line;
        while (std::getline(from, line)) {
            std::string::size_type pos = line.find(pattern);
            if (pos != std::string::npos) {
                std::istringstream iss(line);
                char ignore;
                iss.ignore(std::strlen(pattern), '\n');
                iss >> ignore >> natoms;
                break;
            }
        }
    }
    if (natoms == 0) {
        throw Gauss_error(
            "could not determine number of atoms from Gaussian file");
    }
    from.clear();
    from.seekg(orig_pos, std::ios_base::beg);  // move to original position
    from.clear();

    return natoms;
}

void Gauss_data::get_curr_cart_coord(struct Gauss_coord& coord) const
{
    std::streampos orig_pos = from.tellg();

    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file
    from.clear();

    // Get number of atoms:
    coord.natoms = get_natoms();

    coord.atnum.resize(coord.natoms);
    coord.xyz.resize(coord.natoms, 3);

    std::string line;

    if (filetype == out) {
        throw Gauss_error("not implemented for Gaussian output files");
    }
    else {  // filetype == fchk
        while (std::getline(from, line)) {
            // Get atomic numbers:
            std::string::size_type pos = line.find("Atomic numbers");
            if (pos != std::string::npos) {
                for (int i = 0; i < coord.natoms; ++i) {
                    int atomic_number;
                    from >> atomic_number;
                    coord.atnum(i) = atomic_number;
                }
            }
            // Get current Cartesian coordinates:
            pos = line.find("Current cartesian coordinates");
            if (pos != std::string::npos) {
                for (int i = 0; i < coord.natoms; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        double x;
                        from >> x;
                        coord.xyz(i, j) = x;
                    }
                }
                break;
            }
        }
    }
    from.clear();
    from.seekg(orig_pos, std::ios_base::beg);  // move to original position
    from.clear();
}

void Gauss_data::get_freqs(srs::dvector& freqs) const
{
    if (filetype == fchk) {
        throw Gauss_error("not implemented for Gaussian fchk files");
    }
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    const std::string pattern = " Frequencies --";

    std::string line;
    double v;

    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            std::istringstream iss(line);
            iss.ignore(pattern.size(), '\n');
            while (iss >> v) {
                freqs.push_back(v);
            }
        }
    }
}

void Gauss_data::get_hessians(srs::packed_dmatrix& hess) const
{
    if (filetype == out) {
        throw Gauss_error("not implemented for Gaussian output files");
    }
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    const std::string pattern = "Cartesian Force Constants";

    std::string line;
    std::string buffer;
    int n;

	srs::dvector tmp;
    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            std::istringstream iss(line);
            iss >> buffer >> buffer >> buffer >> buffer >> buffer >> n;
            tmp.resize(n);
            for (int i = 0; i < n; ++i) {
                from >> tmp(i);
                if (!from) {
                    throw Gauss_error("could not read Hessians from fchk file");
                }
            }
        }
    }
    hess = srs::packed_dmatrix(tmp);
}

void Gauss_data::get_pes_scan_data(std::string& scan_coord,
                                   srs::dvector& coord,
                                   srs::dvector& energy) const
{
    if (filetype == fchk) {
        throw Gauss_error("not implemented for Gaussian fchk files");
    }
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    const std::string summary_start
        = "Summary of Optimized Potential Surface Scan";

    scan_coord = get_modredundant_coord();

    std::string line;
    std::string token;
    std::string ignore;
    double value;

    bool summary_found = false;

    while (std::getline(from, line)) {
        if (line.find(summary_start, 0) != std::string::npos) {
            summary_found = true;
            while (std::getline(from, line)) {
                std::istringstream iss(line);
                iss >> token;
                if (token == "Eigenvalues") {
                    iss >> ignore;
                    while (iss >> value) {
                        energy.push_back(value);
                    }
                }
                else if (token == scan_coord) {
                    while (iss >> value) {
                        coord.push_back(value);
                    }
                }
            }
        }
    }
    if (!summary_found) {
        throw Gauss_error(
            "Summary of Optimized Potential Surface Scan not found");
    }
    if (energy.size() != coord.size()) {
        throw Gauss_error("bad number of data read");
    }
}

void Gauss_data::get_nmr_data(std::vector<Gauss_NMR>& nmr,
                              const std::string& nmr_method,
                              const double degen_tol) const
{
    if (filetype == fchk) {
        throw Gauss_error("not implemented for Gaussian fchk files");
    }
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    // Get NMR data:

    std::string line;

    Gauss_NMR_line nmr_line_tmp;
    std::vector<Gauss_NMR_line> nmr_data;

    while (std::getline(from, line)) {
        if (line.find(nmr_method) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find("*****") != std::string::npos) {
                    break;
                }
                int n;
                std::string at;
                std::string buf;
                double sig;

                std::istringstream iss(line);
                if (iss >> n >> at >> buf >> buf >> sig) {
                    nmr_line_tmp.number = n;
                    nmr_line_tmp.atom   = at;
                    nmr_line_tmp.shield = sig;
                    nmr_data.push_back(nmr_line_tmp);
                }
            }
        }
    }
    // Sort shieldings and condense degenerate peaks:

    Gauss_NMR nmr_tmp;

    std::sort(nmr_data.begin(), nmr_data.end());

    for (std::size_t i = 0; i < nmr_data.size(); ++i) {
        nmr_tmp.atom = nmr_data[i].atom;
        nmr_tmp.shield.push_back(nmr_data[i].shield);
        nmr_tmp.number.push_back(nmr_data[i].number);
        if (std::abs(nmr_data[i + 1].shield - nmr_data[i].shield) > degen_tol) {
            nmr.push_back(nmr_tmp);
            nmr_tmp.number.clear();
            nmr_tmp.shield.clear();
        }
    }
}

int Gauss_data::get_no_irc_points() const
{
    std::streampos orig_pos = from.tellg();

    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file
    from.clear();

    int npoints = 0;

    std::string line;
    if (filetype == out) {
        const char pattern[] = "-- Optimized point #";
        while (std::getline(from, line)) {
            if (line.find(pattern, 0) != std::string::npos) {
                std::istringstream iss(srs::trim(line, " "));
                iss.ignore(std::strlen(pattern), '\n');
                iss >> npoints;
            }
        }
        npoints += 1;  // include starting point
    }
    else {  // filetype == fchk
        const char pattern[]
            = "IRC point       1 Results for each geome   R   N=";

        while (std::getline(from, line)) {
            if (line.find(pattern, 0) != std::string::npos) {
                std::istringstream iss(line);
                iss.ignore(std::strlen(pattern), '\n');
                iss >> npoints;
                break;
            }
        }
        npoints /= 2;
    }
    from.clear();
    from.seekg(orig_pos,
               std::ios_base::beg);  // move back to original position
    from.clear();

    return npoints;
}

void Gauss_data::get_irc_data(srs::dvector& mep) const
{
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    std::string line;
    if (filetype == out) {
        const char* pattern;
        if (get_version() == g03) {
            pattern = "Summary of reaction path following:";
        }
        else {
            pattern = "SUMMARY OF REACTION PATH FOLLOWING:";
        }
        int count = 0;
        while (std::getline(from, line)) {
            if (line.find(pattern, 0) != std::string::npos) {
                for (int i = 0; i < 3; ++i) {  // ignore three lines
                    from.ignore(256, '\n');
                }
                const int npoints = get_no_irc_points();
                double vmep;
                double smep;
                int dummy;
                while (count < npoints) {
                    std::getline(from, line);
                    std::istringstream iss(line);
                    iss >> dummy >> vmep >> smep;
                    if (iss) {
                        if (smep != 0.0) {  // do not include starting point
                            mep.push_back(vmep);
                            mep.push_back(smep);
                        }
                        count++;
                    }
                }
            }
        }
        if (count == 0) {
            throw Gauss_error("could not find IRC data in Gaussian file");
        }
    }
    else {  // filetype == fchk
        const char pattern[]
            = "IRC point       1 Results for each geome   R   N=";

        int n = 0;
        while (std::getline(from, line)) {
            if (line.find(pattern, 0) != std::string::npos) {
                std::istringstream iss(line);
                iss.ignore(std::strlen(pattern), '\n');
                iss >> n;
                break;
            }
        }
        if (n > 0) {
            double val;
            for (int i = 0; i < n; ++i) {
                from >> val;
                mep.push_back(val);
            }
        }
        else {
            throw Gauss_error("could not find IRC data in Gaussian file");
        }
    }
}

void Gauss_data::get_irc_geom(srs::dvector& geom) const
{
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    std::string line;
    if (filetype == out) {
        const char pattern_geom[] = "Z-Matrix orientation:";
        const char pattern_opt[]  = "-- Optimized point #";

        const int natoms  = get_natoms();
        const int natoms3 = 3 * natoms;
        const int npoints = get_no_irc_points();

        Gauss_version version = get_version();

        // Temporary store all geometries; save last geometry if an
        // optimized point was found.

        int count = 1;
        srs::dvector geom_tmp(natoms3);
        while (count < npoints) {
            std::getline(from, line);
            if (line.find(pattern_geom, 0) != std::string::npos) {
                geom_tmp.clear();
                for (int i = 0; i < 4; ++i) {  // ignore four lines
                    from.ignore(256, '\n');
                }
                double dum1;
                double dum2;
                double dum3;
                double x;
                double y;
                double z;
                for (int i = 0; i < natoms; ++i) {
                    std::getline(from, line);
                    std::istringstream iss(line);
                    if (version == g94) {
                        iss >> dum1 >> dum2 >> x >> y >> z;
                    }
                    else {
                        iss >> dum1 >> dum2 >> dum3 >> x >> y >> z;
                    }
                    geom_tmp.push_back(x);
                    geom_tmp.push_back(y);
                    geom_tmp.push_back(z);
                }
            }
            if (line.find(pattern_opt, 0) != std::string::npos) {
                for (int i = 0; i < natoms3; ++i) {
                    geom.push_back(geom_tmp(i));
                }
                count++;
            }
        }
        if (count == 0) {
            throw Gauss_error("could not find IRC geometries in Gaussian file");
        }
    }
    else {  // filetype == fchk
        const char pattern[]
            = "IRC point       1 Geometries               R   N=";

        int n = 0;
        while (std::getline(from, line)) {
            if (line.find(pattern, 0) != std::string::npos) {
                std::istringstream iss(line);
                iss.ignore(std::strlen(pattern), '\n');
                iss >> n;
                break;
            }
        }
        if (n > 0) {
            double val;
            for (int i = 0; i < n; ++i) {
                from >> val;
                geom.push_back(val);
            }
        }
        else {
            throw Gauss_error("could not find IRC geometries in Gaussian file");
        }
    }
}

void Gauss_data::get_irc_grad(srs::dvector& grad) const
{
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    std::string line;
    if (filetype == out) {
        const char pattern_grad[]
            = "Center     Atomic                   Forces";
        const char pattern_opt[] = "-- Optimized point #";

        const int natoms  = get_natoms();
        const int natoms3 = 3 * natoms;
        const int npoints = get_no_irc_points();

        // Temporary store all gradients; save last gradient if an
        // optimized point was found.

        int count = 1;
        srs::dvector grad_tmp(natoms3);
        while (count < npoints) {
            std::getline(from, line);
            if (line.find(pattern_grad, 0) != std::string::npos) {
                grad_tmp.clear();
                for (int i = 0; i < 2; ++i) {  // ignore two lines
                    from.ignore(256, '\n');
                }
                double dum1;
                double dum2;
                double dx;
                double dy;
                double dz;
                for (int i = 0; i < natoms; ++i) {
                    std::getline(from, line);
                    std::istringstream iss(line);
                    iss >> dum1 >> dum2 >> dx >> dy >> dz;
                    grad_tmp.push_back(dx);
                    grad_tmp.push_back(dy);
                    grad_tmp.push_back(dz);
                }
            }
            if (line.find(pattern_opt, 0) != std::string::npos) {
                for (int i = 0; i < natoms3; ++i) {
                    grad.push_back(-1.0 * grad_tmp(i));  // forces to gradients
                }
                count++;
            }
        }
        if (count == 0) {
            throw Gauss_error("could not find IRC gradients in Gaussian file");
        }
    }
    else {  // filetype == fchk
        const char pattern[]
            = "IRC point       1 Gradient at each geome   R   N=";

        int n = 0;
        while (std::getline(from, line)) {
            if (line.find(pattern, 0) != std::string::npos) {
                std::istringstream iss(line);
                iss.ignore(std::strlen(pattern), '\n');
                iss >> n;
                break;
            }
        }
        if (n > 0) {
            double val;
            for (int i = 0; i < n; ++i) {
                from >> val;
                grad.push_back(val);
            }
        }
        else {
            throw Gauss_error("could not find IRC gradients in Gaussian file");
        }
    }
}

void Gauss_data::get_irc_hess(srs::dvector& hess) const
{
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    if (filetype == out) {
        const char pattern_hess[] = " The second derivative matrix:";
        const char pattern_opt[]  = "-- Optimized point #";

        const int npoints = get_no_irc_points();
        const int natoms  = get_natoms();
        const int natoms3 = 3 * natoms;
        const int nhess   = natoms3 * (natoms3 + 1) / 2;

        // Temporary store all Hessians; save last Hessians if an
        // optimized point was found.

        int count = 1;
        double fc;
        std::string line;
        std::string token;
        std::string data;
        srs::dvector hess_tmp(nhess);
        while (count < npoints) {
            std::getline(from, line);
            if (line.find(pattern_hess, 0) != std::string::npos) {
                hess_tmp.clear();
                while (from >> token) {
                    std::getline(from, data);
                    std::istringstream iss(data);
                    while (iss >> fc) {
                        hess_tmp.push_back(fc);
                    }
                    if (hess_tmp.size() >= nhess) {
                        break;
                    }
                }
            }
            if (line.find(pattern_opt, 0) != std::string::npos) {
                for (int i = 0; i < nhess; ++i) {
                    hess.push_back(hess_tmp(i));
                }
                count++;
            }
        }
        if (count == 0) {
            throw Gauss_error("could not find IRC Hessians in Gaussian file");
        }
    }
    else {  // filetype == fchk
        throw Gauss_error(
            "IRC Hessians can only be extracted from Gaussian output "
            "files");
    }
}

std::string Gauss_data::get_modredundant_coord() const
{
    if (filetype == fchk) {
        throw Gauss_error("not implemented for Gaussian fchk files");
    }
    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file

    const std::string pattern = "!    Initial Parameters    !";

    std::string line;
    std::string name;
    std::string def;
    std::string deriv;
    double val;
    char ch;

    while (std::getline(from, line)) {
        if (line.find(pattern, 0) != std::string::npos) {
            for (int i = 0; i < 4; ++i) {
                std::getline(from, line);  // ignore
            }
            while (std::getline(from, line)) {
                std::istringstream iss(line);
                iss >> ch >> name >> def >> val >> deriv;
                if (deriv == "Scan") {
                    return name;
                }
            }
        }
    }
    throw Gauss_error("ModRedundant coordinate not found");
}
