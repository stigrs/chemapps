#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdexcept>
#include <exception>

#include <chem/input.h>
#include <chem/utils.h>
#include <armadillo>


int main(int /* argc */, char* argv[])
{
    int i;
    double d;
    std::string s;
    
    arma::ivec iv;
    arma::uvec uv;
    arma::vec  dv;
    arma::ivec iv_ans = { 1, 2, 3, 4 };
    arma::uvec uv_ans = { 5, 6, 7, 8 };
    arma::vec  dv_ans = { 0.1, 0.2, 0.3, 0.4, 0.5 };

    std::map<std::string, Input> data;
    typedef std::map<std::string, Input>::iterator Input_iter;

    data["integer"] = Input(i);
    data["double"]  = Input(d);
    data["string"]  = Input(s);
    data["ivector"] = Input(iv);
    data["uvector"] = Input(uv);
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
        for (arma::uword it = 0; it < iv_ans.size(); ++it) {
            chem::Assert(iv(it) == iv_ans(it), 
                         std::runtime_error("ivector failed"));
        }
        for (arma::uword it = 0; it < uv_ans.size(); ++it) {
            chem::Assert(uv(it) == uv_ans(it), 
                         std::runtime_error("uvector failed"));
        }
        for (arma::uword it = 0; it < dv_ans.size(); ++it) {
            chem::Assert(dv(it) == dv_ans(it), 
                         std::runtime_error("dvector failed"));
        }
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
