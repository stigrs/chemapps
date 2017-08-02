#include <stdexcept>
#include <exception>

#include <chem/datum.h>
#include <chem/utils.h>

int main(int /* argc */, char* argv[])
{
    try {
        chem::Assert(datum::h == 6.62607004000e-34,
                     std::runtime_error("wrong Planck constant"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
