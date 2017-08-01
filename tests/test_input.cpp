#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdexcept>
#include <exception>

#include <chem/input.h>
#include <chem/utils.h>
#include <armadillo>


int main(int argc, char* argv[])
{
    int i;
    double d;
    std::string s;
    
    arma::ivec iv;
    arma::vec  dv;
    arma::ivec iv_ans = { 1, 2, 3, 4 };
    arma::vec  dv_ans = { 0.1, 0.2, 0.3, 0.4, 0.5 };

    std::map<std::string, Input> data;
    typedef std::map<std::string, Input>::iterator Input_iter;

    data["integer"] = Input(i);
    data["double"]  = Input(d);
    data["string"]  = Input(s);
    data["ivector"] = Input(iv);
    data["dvector"] = Input(dv);

    std::ifstream from;

    try {
        chem::fopen(from, "test_input.inp");

        std::string key;
        while (from >> key) {
            Input_iter it = data.find(key);
            if (it != data.end()) {
                from >> it->second;
            }
        }

        chem::Assert(i == 1, std::runtime_error("integer failed"));
        chem::Assert(d == 2.0, std::runtime_error("double failed"));
        chem::Assert(s == "hello", std::runtime_error("string failed"));
        for (int i = 0; i < iv_ans.size(); ++i) {
            chem::Assert(iv(i) == iv_ans(i), 
                         std::runtime_error("ivector failed"));
        }
        for (int i = 0; i < dv_ans.size(); ++i) {
            chem::Assert(dv(i) == dv_ans(i), 
                         std::runtime_error("dvector failed"));
        }
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
