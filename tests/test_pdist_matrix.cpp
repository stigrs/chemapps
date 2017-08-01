#include <iostream>
#include <stdexcept>
#include <exception>
#include <chem/math.h>
#include <chem/utils.h>
#include <armadillo>


int main(int argc, char* argv[])
{
    try {
        arma::mat mat(4,3);
        mat(0,0) = 0.0;
        mat(0,1) = 0.0;
        mat(0,2) = 0.0;
        mat(1,0) = 1.0;
        mat(1,1) = 1.0;
        mat(1,2) = 1.0;
        mat(2,0) = 2.0;
        mat(2,1) = 2.0;
        mat(2,2) = 2.0;
        mat(3,0) = 3.0;
        mat(3,1) = 3.0;
        mat(3,2) = 3.0;
        
        arma::mat dm;
        chem::pdist_matrix(dm, mat);
        
        arma::mat dm_ans = arma::zeros<arma::mat>(4,4);
        
        dm_ans(0,1) = 1.73205081;
        dm_ans(0,2) = 3.46410162;
        dm_ans(0,3) = 5.19615242;
        dm_ans(1,2) = 1.73205081;
        dm_ans(1,3) = 3.46410162;
        dm_ans(2,3) = 1.73205081;
        dm_ans(1,0) = dm_ans(0,1);
        dm_ans(2,0) = dm_ans(0,2);
        dm_ans(3,0) = dm_ans(0,3);
        dm_ans(2,1) = dm_ans(1,2);
        dm_ans(3,1) = dm_ans(1,3);
        dm_ans(3,2) = dm_ans(2,3);

        for (int j = 0; j < dm.n_cols; ++j) {
            for (int i = j; i < dm.n_rows; ++i) {
                if (i != j) {
                    chem::Assert(std::abs(dm(i,j) - dm_ans(i,j)) < 1.0e-8, 
                                 std::runtime_error("bad pdist_matrix"));
                }
            }
        }
    } 
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
    }
}
