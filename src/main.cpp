#include <iostream>
#include "namu/tree_summary.hpp"

using namespace namu;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    std::cout << "Starting..." << std::endl;

    TreeSummary sumt;

    try {
        sumt.read_tree_file("../data/test.tre", 1);
    }
    catch(const std::exception & x) {
        std::cerr << "Program aborting due to errors encountered reading tree file." << std::endl;
        std::cerr << x.what() << std::endl;
        std::exit(0);
    }

    sumt.show_summary();

    std::cout << "\nFinished!" << std::endl;

    return 0;
}
