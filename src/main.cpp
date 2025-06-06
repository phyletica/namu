#include <iostream>
#include "namu/node.hpp"
#include "namu/tree.hpp"
#include "namu/tree_manip.hpp"
#include "namu/tree_io.hpp"

using namespace namu;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    std::cout << "Starting..." << std::endl;
    TreeManip::SharedPtrVector trees;
    /* std::string newick = "(1:0.3,2:0.3,(3:0.2,(4:0.1,5:0.1):0.1):0.1);"; */
    std::string newick = "(1[&height=0]:0.3,2[&height=0]:0.3,(3[&height=0]:0.2,(4[&height=0]:0.1,5[&height=0]:0.1)[&height=0.1,height_index=0]:0.1)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2];";
    std::cout << "Input: " << newick << std::endl;
    TreeIO tree_io = TreeIO();
    trees = tree_io.parse_from_newick(newick);
    std::cout << "Number of trees parsed: " << trees.size() << std::endl;
    TreeManip::SharedPtr tm = trees.at(0);
    std::cout << "Output: " << tm->make_newick(3) << std::endl;
    std::cout << "\nFinished!" << std::endl;

    return 0;
}
