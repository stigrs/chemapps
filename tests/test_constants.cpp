#include <stdexcept>
#include <exception>

#include <chem/constants.h>
#include <chem/utils.h>

int main(int argc, char* argv[])
{
    try {
        double h = constants::planck;
        chem::Assert(h == 6.62607004000e-34,
                     std::runtime_error("wrong Planck constant"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
