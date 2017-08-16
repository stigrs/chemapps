/**
   @file mopac.cpp

   This file is part of ChemApps - A C++ Chemistry Toolkit

   Copyright (C) 2016-2017  Stig Rune Sellevag

   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <chem/datum.h>
#include <chem/input.h>
#include <chem/mopac.h>
#include <chem/utils.h>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

Mopac::Mopac(std::istream& from, const std::string& key)
{
    typedef std::map<std::string, Input>::iterator Input_iter;
    typedef std::map<std::string, Input>::const_iterator Cinput_iter;

    // Read input data:

    std::string version_def  = "mopac5022mn";
    std::string keywords_def = "PM6-D GEO-OK PRECISE";
    std::string jobname_def  = "mopac";
    int opt_geom_def         = 1;

    std::map<std::string, Input> input_data;
    input_data["version"]  = Input(version, version_def);
    input_data["jobname"]  = Input(jobname, jobname_def);
    input_data["opt_geom"] = Input(opt_geom, opt_geom_def);

    if (chem::find_section(from, key)) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "keywords") {
                std::string line;
                std::getline(from, line);
                if (line.empty()) {  // not entirely safe
                    std::getline(from, line);
                }
                keywords = chem::trim(line, " ");
            }
            else {
                Input_iter it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }
    else {
        throw Mopac_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (Cinput_iter it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Mopac_error(it->first + " not initialized");
        }
    }
}

void Mopac::run(Molecule& mol) const
{
    // Create Mopac input file:
    write_dat(mol);

    // Run Mopac:
    std::string cmd = version + " " + jobname + ".dat";
    if (std::system(cmd.c_str()) != 0) {
        throw Mopac_error("running " + version + " failed");
    }

    // Check convergence:
    if (!check_convergence()) {
        throw Mopac_error(jobname + " failed to converge");
    }

    // Update molecular energy:
    mol.set_elec_energy(get_heat_of_formation());

    // Update Cartesian coordinates:
    arma::mat xyz = mol.get_xyz();
    get_xyz(xyz);
    mol.set_xyz(xyz);
}

void Mopac::write_dat(const Molecule& mol) const
{
    std::ofstream to;
    chem::fopen(to, jobname + ".dat");
    to << keywords << '\n' << mol.get_title() << "\n\n";
    write_xyz(to, mol);
}

void Mopac::write_xyz(std::ostream& to, const Molecule& mol) const
{
    chem::Format<double> fix;
    fix.fixed().width(10).precision(6);

    for (std::size_t i = 0; i < mol.get_atoms().size(); ++i) {
        to << mol.get_atoms()[i].atomic_symbol << '\t';
        for (arma::uword j = 0; j < mol.get_xyz().n_cols; ++j) {
            to << fix(mol.get_xyz()(i, j)) << " " << opt_geom << " ";
        }
        to << '\n';
    }
}

bool Mopac::check_convergence() const
{
    bool converged = false;

    std::ifstream from;
    chem::fopen(from, jobname + ".out");

    std::string buf;
    while (std::getline(from, buf)) {
        std::string pattern = "SCF FIELD WAS ACHIEVED";
        if (buf.find(pattern) != std::string::npos) {
            converged = true;
            break;
        }
    }
    return converged;
}

double Mopac::get_heat_of_formation() const
{
    double heat = 0.0;
    bool found  = false;

    std::ifstream from;
    chem::fopen(from, jobname + ".out");

    std::string buf;
    while (std::getline(from, buf)) {
        std::string pattern = "SCF FIELD WAS ACHIEVED";  // is this fail-safe?
        if (buf.find(pattern) != std::string::npos) {
            while (std::getline(from, buf)) {
                pattern = "FINAL HEAT OF FORMATION = ";
                if (buf.find(pattern) != std::string::npos) {
                    std::istringstream iss(buf);
                    for (int i = 0; i < 5; ++i) {  // ignore words in pattern
                        iss >> buf;
                    }
                    iss >> heat;
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
    }
    if (!found) {
        throw Mopac_error("final heat of formation not found");
    }
    else {
        return heat * datum::cal_to_J;
    }
}

void Mopac::get_xyz(arma::mat& xyz) const
{
    // Note: xyz must have the correct size on input, no resizing is done.

    double heat = 0.0;
    bool found  = false;

    std::ifstream from;
    chem::fopen(from, jobname + ".out");

    std::string buf;
    while (std::getline(from, buf)) {
        std::string pattern = "SCF FIELD WAS ACHIEVED";  // is this fail-safe?
        if (buf.find(pattern) != std::string::npos) {
            while (std::getline(from, buf)) {
                pattern = "CARTESIAN COORDINATES";
                if (buf.find(pattern) != std::string::npos) {
                    for (int i = 0; i < 3; ++i) {  // ignore three lines
                        std::getline(from, buf);
                    }
                    int center;
                    std::string atom;
                    double x;
                    double y;
                    double z;
                    for (arma::uword i = 0; i < xyz.n_rows; ++i) {
                        std::getline(from, buf);
                        std::istringstream iss(buf);
                        iss >> center >> atom >> x >> y >> z;
                        xyz(i, 0) = x;
                        xyz(i, 1) = y;
                        xyz(i, 2) = z;
                    }
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
    }
    if (!found) {
        throw Mopac_error("optimized Cartesian coordinates not found");
    }
}
