#include <iostream>
#include <fstream>
#include <stdexcept>
#include <exception>

#include <chem/element.h>
#include <chem/utils.h>
#include <chem/molecule.h>
#include <armadillo>


int main(int /* argc */, char* argv[])
{
    arma::mat xyz_ans(8,3);
    xyz_ans(0,0) =  0.0000;
    xyz_ans(0,1) =  0.0000;
    xyz_ans(0,2) =  0.7637;
    xyz_ans(1,0) =  0.0000;
    xyz_ans(1,1) =  0.0000;
    xyz_ans(1,2) = -0.7637;
    xyz_ans(2,0) =  0.0000;
    xyz_ans(2,1) =  1.0121;
    xyz_ans(2,2) =  1.1564;
    xyz_ans(3,0) = -0.8765;
    xyz_ans(3,1) = -0.5060;
    xyz_ans(3,2) =  1.1564;
    xyz_ans(4,0) =  0.8765;
    xyz_ans(4,1) = -0.5060;
    xyz_ans(4,2) =  1.1564;
    xyz_ans(5,0) =  0.0000;
    xyz_ans(5,1) = -1.0121;
    xyz_ans(5,2) = -1.1564;
    xyz_ans(6,0) = -0.8765;
    xyz_ans(6,1) =  0.5060;
    xyz_ans(6,2) = -1.1564;
    xyz_ans(7,0) =  0.8765;
    xyz_ans(7,1) =  0.5060;
    xyz_ans(7,2) = -1.1564;


    try {
        std::ifstream from;
        chem::fopen(from, "test_molecule.inp");

        Molecule mol(from, std::cout, "Molecule");
        arma::mat xyz = mol.get_xyz();
        
        chem::Assert(arma::approx_equal(xyz, xyz_ans, "absdiff", 1.0e-12), 
                     std::runtime_error("bad xyz coords"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
