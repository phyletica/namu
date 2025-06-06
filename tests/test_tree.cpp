#include <catch2/catch.hpp>
#include "namu/tree.hpp"
#include "namu/tree_manip.hpp"
#include "namu/tree_io.hpp"

const double namu::Node::_smallest_edge_length = 1.0e-12;

TEST_CASE("Testing Tree bare constructor", "[Tree]") {
    namu::Tree t = namu::Tree();
    REQUIRE(t.num_nodes() == 0);
}

TEST_CASE("Testing strom tutorial newick", "[Tree]") {
    namu::TreeManip::SharedPtrVector trees;
    std::string newick = "(1:0.3,2:0.3,(3:0.2,(4:0.1,5:0.1):0.1):0.1);";
    std::cout << "Input: " << newick << std::endl;
    namu::TreeIO tree_io = namu::TreeIO();
    trees = tree_io.parse_from_newick(newick);
    std::cout << "Number of trees parsed: " << trees.size() << std::endl;
    REQUIRE(trees.size() == 1);
    namu::TreeManip::SharedPtr tm = trees.at(0);
    std::string new_str = tm->make_newick(3);
    std::cout << "Output: " << new_str << std::endl;
    std::string expected_new_str = "(1:0.300,2:0.300,(3:0.200,(4:0.100,5:0.100):0.100):0.100)";
    REQUIRE(new_str == expected_new_str);
}
