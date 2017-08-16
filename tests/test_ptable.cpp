#include <chem/ptable.h>
#include <chem/utils.h>
#include <exception>
#include <stdexcept>

int main(int /* argc */, char* argv[])
{
    try {
        double atomic_weight = ptable::get_element("C").atomic_weight;
        chem::Assert(atomic_weight == (12.0096 + 12.0116) * 0.5,
                     std::runtime_error("wrong atomic weight for C"));

        double atomic_mass = ptable::get_element("109Ag").atomic_mass;
        chem::Assert(atomic_mass == 108.904755,
                     std::runtime_error("wrong atomic mass for 109Ag"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
