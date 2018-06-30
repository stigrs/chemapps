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
#include <chem/gauss_error.h>
#include <chem/gaussian.h>
#include <srs/utils.h>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl>
#include <map>

void Gaussian::init(std::istream& from, const std::string& key)
{
    // Read input data:

    int nosave_tmp;
    std::map<std::string, srs::Input> input_data;
    input_data["version"]     = srs::Input(version, "rung09");
    input_data["jobname"]     = srs::Input(jobname, "gauss");
    input_data["nprocshared"] = srs::Input(nprocshared, 1);
    input_data["nosave"]      = srs::Input(nosave_tmp, 1);

    if (srs::find_section(from, key)) {
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
                keywords = srs::trim(line, " ");
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Gauss_error(it->first + " not initialized");
        }
    }

    // Validate input data:

    Expects(nosave_tmp == 0 || nosave_tmp == 1);
    nosave = false;
    if (nosave_tmp == 1) {
        nosave = true;
    }
    Expects(nprocshared >= 1);
}

void Gaussian::run(Molecule& mol) const
{
    // Create Gaussian input file:
    write_com(mol);

    // Run Gaussian:
    std::string cmd = version + " " + jobname;
    if (std::system(cmd.c_str()) != 0) {
        throw Gauss_error("running " + version + " failed");
    }

    // Get Gaussian output data:
    std::ifstream logfile;
    srs::fopen(logfile, jobname + ".log");

    Gauss_data data(logfile, out);

    // Check termination:
    if (!data.check_termination()) {
        throw Gauss_error("Gaussian did not terminate normally");
    }

    // Check geometry convergence:
    if (!data.check_opt_conv()) {
        throw Gauss_error("stationary point not found");
    }

    // Update molecular energy:
    mol.set_elec_energy(data.get_scf_zpe_energy()[0]);

    // Update Cartesian coordinates:
    Gauss_coord coord;
    data.get_opt_cart_coord(coord);
    mol.set_xyz(coord.xyz);
}

//------------------------------------------------------------------------------

void Gaussian::write_com(const Molecule& mol) const
{
    std::ofstream to;
    srs::fopen(to, jobname + ".com");
    to << "%nprocshared=" << nprocshared << '\n'
       << "%chk=" << jobname << ".chk" << '\n';
    if (nosave) {
        to << "%nosave\n";
    }
    to << "# " << keywords << "\n\n"
       << mol.get_title() << "\n\n"
       << mol.get_charge() << " " << mol.get_elec_state()[0] << '\n';

    srs::Format<double> fix;
    fix.fixed().width(10).precision(6);
    for (std::size_t i = 0; i < mol.get_atoms().size(); ++i) {
        to << mol.get_atoms()[i].atomic_symbol << '\t';
        for (srs::size_t j = 0; j < mol.get_xyz().cols(); ++j) {
            to << fix(mol.get_xyz()(i, j)) << "  ";
        }
        to << '\n';
    }
    to << '\n';
}
