#include <iostream>
#include "namu/namu.hpp"

using namespace namu;

// static data member initializations
std::string  Namu::_program_name        = "namu";
unsigned     Namu::_major_version       = 0;
unsigned     Namu::_minor_version       = 1;
const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    Namu namu;

    try {
        namu.process_command_line_options(argc, argv);
        namu.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception caught: " << std::endl << typeid(x).name() << std::endl;
        std::cerr << "Exception message: " << std::endl << x.what() << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
}
