// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/periodic_table.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

// Error reporting:

struct Setup_error : std::runtime_error {
    Setup_error(const std::string& s) : std::runtime_error(s) {}
};

struct Init_geom_error : std::runtime_error {
    Init_geom_error(const std::string& s) : std::runtime_error(s) {}
};

struct Opt_geom_error : std::runtime_error {
    Opt_geom_error(const std::string& s) : std::runtime_error(s) {}
};

struct Create_input_error : std::runtime_error {
    Create_input_error(const std::string& s) : std::runtime_error(s) {}
};

struct Orbitals_error : std::runtime_error {
    Orbitals_error(const std::string& s) : std::runtime_error(s) {}
};

struct Gms_run_error : std::runtime_error {
    Gms_run_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Declarations:

enum Runtype_t { energy, optimize, sadpoint };
enum Orbital_t { guess, mcscf_opt };
enum Coord_t { cart, zmat };
enum Axis_t { x_axis, y_axis, z_axis };

void parse_setup(std::ifstream& from);
void init_geom();
void get_opt_geom(const std::string& filename);
void print_geom(std::ostream& to);
void step_geom();
void create_input(const std::string& new_inp_file,
                  const std::string& old_dat_file);
void get_orbitals(const std::string& filename);
void summarize();
bool gms_normal_exit(const std::string& filename);
double gms_final_energy(const std::string& filename);

//-----------------------------------------------------------------------------

// Global variables:

Numlib::Vec<int> fragment;
std::vector<double> at_num;
std::vector<std::string> at_sym;
std::vector<std::string> orbitals;
std::vector<std::string> log_files;
Numlib::Mat<double> geom;

std::string progname;
std::string tmlfile;
std::string jobname;
Runtype_t runtype;
Orbital_t orbstart;
Coord_t coord;
Axis_t axis;
double stepsize;
int nstep;

//-----------------------------------------------------------------------------

// Program for performing surface scan calculations with GAMESS.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " setup_file\n";
        return 1;
    }

    try {
        std::ifstream from;
        Stdutils::fopen(from, args[1]);

        parse_setup(from);
        init_geom();

        // Perform GAMESS surface scan:

        int status;
        std::string cmd;
        std::string inp_file;
        std::string dat_file = tmlfile;

        for (int i = 0; i < nstep; i++) {
            std::ostringstream oss;
            oss << jobname << "_" << i + 1;
            inp_file = oss.str() + ".inp";
            log_files.push_back(oss.str() + ".log");

            create_input(inp_file.c_str(), dat_file.c_str());

            cmd = progname + " " + inp_file;
            status = std::system(cmd.c_str()); // run GAMESS
            if (status != 0) {
                throw Gms_run_error("running " + cmd + " failed");
            }
            if (!gms_normal_exit(log_files[i].c_str())) {
                throw Gms_run_error("execution of GAMESS failed, check " +
                                    log_files[i]);
            }
            if ((runtype == optimize) && (coord == cart)) {
                get_opt_geom(log_files[i].c_str());
            }
            dat_file = oss.str() + ".dat"; // use orbitals from previous calc.
        }

        // Summarize results:

        summarize();
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Initialize variables by parsing setup file.
void parse_setup(std::ifstream& from)
{
    using namespace Stdutils;

    std::string runtype_tmp;
    std::string orbstart_tmp;
    std::string coord_tmp;
    std::string axis_tmp;
    fragment = {0};

    get_token_value(from, 0, "progname", progname, std::string("gmss"));
    get_token_value(from, 0, "tmlfile", tmlfile);
    get_token_value(from, 0, "runtype", runtype_tmp);
    get_token_value(from, 0, "jobname", jobname);
    get_token_value(from, 0, "coord", coord_tmp);
    get_token_value(from, 0, "orbstart", orbstart_tmp, std::string("guess"));
    get_token_value(from, 0, "axis", axis_tmp, std::string("x_axis"));
    get_token_value(from, 0, "fragment", fragment);
    get_token_value(from, 0, "nstep", nstep);
    get_token_value(from, 0, "stepsize", stepsize);

    // Check if data are sensible:

    if (runtype_tmp == "energy") {
        runtype = energy;
    }
    else if (runtype_tmp == "optimize") {
        runtype = optimize;
    }
    else if (runtype_tmp == "sadpoint") {
        runtype = sadpoint;
    }
    else {
        throw Setup_error("runtype has bad value: " + runtype_tmp);
    }

    if (orbstart_tmp == "guess") {
        orbstart = guess;
    }
    else if (orbstart_tmp == "mcscf_opt") {
        orbstart = mcscf_opt;
    }
    else {
        throw Setup_error("orbstart has bad value: " + orbstart_tmp);
    }

    if (coord_tmp == "cart") {
        coord = cart;
    }
    else if (coord_tmp == "zmat") {
        coord = zmat;
    }
    else {
        throw Setup_error("coord has bad value: " + coord_tmp);
    }

    if (axis_tmp == "x_axis") {
        axis = x_axis;
    }
    else if (axis_tmp == "y_axis") {
        axis = y_axis;
    }
    else if (axis_tmp == "z_axis") {
        axis = z_axis;
    }
    else {
        throw Setup_error("axis has bad value: " + axis_tmp);
    }

    if (nstep < 1) {
        throw Setup_error("nstep has bad value: " + std::to_string(nstep));
    }
}

// Get starting geometry.
void init_geom()
{
    namespace Pt = Chem::Periodic_table;

    std::ifstream from;
    Stdutils::fopen(from, tmlfile.c_str());

    std::string key = "Geometry";

    auto pos = Stdutils::find_token(from, key);
    if (pos != -1) {
        std::string token;
        std::vector<double> geom_tmp;

        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (Pt::atomic_symbol_is_valid(token)) {
                double znuc, x, y, z;
                from >> znuc >> x >> y >> z;
                if (!from) {
                    throw Init_geom_error("cannot read starting geometry");
                }
                geom_tmp.push_back(x);
                geom_tmp.push_back(y);
                geom_tmp.push_back(z);
                at_sym.push_back(token);
                at_num.push_back(Pt::get_atomic_number(token));
            }
        }

        std::size_t natoms = at_num.size();
        if (natoms < 1) {
            throw Init_geom_error("geometry not initialized");
        }

        geom.resize(natoms, 3);
        for (std::size_t i = 0; i < natoms; i++) {
            for (int j = 0; j < 3; j++) {
                geom(i, j) = geom_tmp[i * 3 + j];
            }
        }
    }
    else {
        throw Init_geom_error("cannot find " + key + " section");
    }
}

// Get optimized geometry.
void get_opt_geom(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    std::string key1;
    switch (runtype) {
    case optimize:
        key1 = "***** EQUILIBRIUM GEOMETRY LOCATED *****";
        break;
    case sadpoint:
        key1 = "***** SADDLE POINT LOCATED *****";
        break;
    default:
        break;
    }
    std::string key2 = "COORDINATES OF ALL ATOMS ARE";

    bool found = false;

    std::string line, name;
    double znuc, x, y, z;
    std::vector<double> geom_tmp;

    while (std::getline(from, line)) {
        std::string::size_type pos = line.find(key1);
        if (pos != std::string::npos) { // geometry found
            while (std::getline(from, line)) {
                pos = line.find(key2);
                if (pos != std::string::npos) { // coord. of all atoms found
                    found = true;
                    std::getline(from, line); // ignore two lines
                    std::getline(from, line);
                    for (std::size_t i = 0; i < at_num.size(); i++) {
                        std::getline(from, line);
                        std::istringstream iss(line);
                        iss >> name >> znuc >> x >> y >> z;
                        if (!iss) {
                            throw Opt_geom_error(
                                "cannot read optimized geometry");
                        }
                        geom_tmp.push_back(x);
                        geom_tmp.push_back(y);
                        geom_tmp.push_back(z);
                    }
                }
            }
        }
    }
    if (!found) {
        throw Opt_geom_error("optimized geometry not found in " + filename);
    }

    for (std::size_t i = 0; i < at_num.size(); i++) {
        for (int j = 0; j < 3; j++) {
            geom(i, j) = geom_tmp[i * 3 + j];
        }
    }
}

// Print geometry to output stream.
void print_geom(std::ostream& to)
{
    Stdutils::Format<double> fix1;
    Stdutils::Format<double> fix10;
    fix1.fixed().width(3).precision(1);
    fix10.fixed().width(15).precision(10);

    for (std::size_t i = 0; i < at_num.size(); i++) {
        to << at_sym[i] << "   " << fix1(at_num[i]) << "   ";
        for (int j = 0; j < 3; j++) {
            to << fix10(geom(i, j)) << "   ";
        }
        to << '\n';
    }
}

// Generate new geometry.
void step_geom()
{
    switch (coord) {
    case cart:
        for (auto i : fragment) {
            geom(i - 1, axis) += stepsize;
        }
        break;
    case zmat:
        break;
    }
}

// Create new GAMESS input file.
void create_input(const std::string& new_inp_file,
                  const std::string& old_dat_file)
{
    std::ifstream from;
    std::ofstream to;

    Stdutils::fopen(from, tmlfile);
    Stdutils::fopen(to, new_inp_file);

    std::string key = "Input";

    auto pos = Stdutils::find_token(from, key);
    if (pos != -1) {
        std::string line;
        while (std::getline(from, line)) {
            if (Stdutils::trim(line, " \t") == "End") {
                break;
            }
            else if (Stdutils::trim(line, " \t") == "GEOMETRY_HERE") {
                step_geom();
                print_geom(to);
            }
            else {
                to << line << '\n';
            }
        }
    }
    else {
        throw Create_input_error("cannot find " + key + " section");
    }

    switch (orbstart) {
    case mcscf_opt:
        get_orbitals(old_dat_file);
        for (const auto& orb : orbitals) {
            to << orb << '\n';
        }
    case guess:
    default:
        return;
    }
}

// Get molecular orbitals from previous calculation.
void get_orbitals(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    std::string pattern;
    switch (orbstart) {
    case mcscf_opt:
        pattern = "--- OPTIMIZED MCSCF MO-S ---";
        break;
    default:
        return;
    }

    orbitals.clear();
    std::string line;
    while (std::getline(from, line)) {
        std::string::size_type pos = line.find(pattern);
        if (pos != std::string::npos) { // get final orbitals
            orbitals.clear();
            orbitals.push_back(line);
            while (std::getline(from, line)) {
                orbitals.push_back(line);
                if (Stdutils::trim(line, " \t") == "$END") {
                    break;
                }
            }
        }
    }
    if (orbitals.empty()) {
        throw Orbitals_error("no orbitals found in " + filename);
    }
}

void summarize()
{
    Stdutils::Format<char> line;
    Stdutils::Format<double> fix6;
    Stdutils::Format<double> fix10;
    line.width(40).fill('-');
    fix6.fixed();
    fix10.fixed().width(16).precision(10);

    double s = stepsize;

    std::cout << "No.\ts/angstrom\tEnergy/hartree\n" << line('-') << '\n';
    for (int i = 0; i < nstep; i++) {
        std::cout << i + 1 << '\t' << fix6(s) << '\t'
                  << fix10(gms_final_energy(log_files[i].c_str())) << '\n';
        s += stepsize;
    }
}

bool gms_normal_exit(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    const std::string pattern = "EXECUTION OF GAMESS TERMINATED NORMALLY";

    std::string line;
    while (std::getline(from, line)) {
        std::string::size_type pos = line.find(pattern);
        if (pos != std::string::npos) {
            return true;
        }
    }
    return false;
}

double gms_final_energy(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    const std::string pattern = "TOTAL ENERGY =";

    std::string line, str;
    double energy = 0.0;
    while (std::getline(from, line)) {
        std::string::size_type pos = line.find(pattern);
        if (pos != std::string::npos) {
            str = line.substr(pos + pattern.size());
            std::istringstream iss(str);
            iss >> energy;
        }
    }
    return energy;
}
