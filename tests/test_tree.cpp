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

    namu::Tree::SharedPtr tree = tm->get_tree();
    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);

    REQUIRE(tm->calc_tree_length() == Approx(1.2));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    tm->scale_all_internal_heights(2.0);
    REQUIRE(tm->calc_tree_length() == Approx(2.4));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);
}

TEST_CASE("Testing strom tutorial newick with height comments", "[Tree]") {
    namu::TreeManip::SharedPtrVector trees;
    std::string newick = "(1[&height=0]:0.3,2[&height=0]:0.3,(3[&height=0]:0.2,(4[&height=0]:0.1,5[&height=0]:0.1)[&height=0.1,height_index=0]:0.1)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2];";
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

    namu::Tree::SharedPtr tree = tm->get_tree();
    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);

    REQUIRE(tm->calc_tree_length() == Approx(1.2));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    tm->scale_all_internal_heights(2.0);
    REQUIRE(tm->calc_tree_length() == Approx(2.4));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);
}

TEST_CASE("Testing old tips", "[Tree]") {
    namu::TreeManip::SharedPtrVector trees;
    std::string newick = "(1:0.3,2:0.25,(3:0.2,(4:0.08,5:0.1):0.1):0.1);";
    namu::TreeIO tree_io = namu::TreeIO();
    trees = tree_io.parse_from_newick(newick);
    REQUIRE(trees.size() == 1);
    namu::TreeManip::SharedPtr tm = trees.at(0);
    std::string new_str = tm->make_newick(3);
    std::string expected_new_str = "(1:0.300,2:0.250,(3:0.200,(4:0.080,5:0.100):0.100):0.100)";
    REQUIRE(new_str == expected_new_str);

    namu::Tree::SharedPtr tree = tm->get_tree();
    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);

    // 1.2 - 0.07
    REQUIRE(tm->calc_tree_length() == Approx(1.13));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    tm->scale_all_internal_heights(2.0);
    // 2.4 - 0.07
    REQUIRE(tm->calc_tree_length() == Approx(2.33));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);
}

TEST_CASE("Testing old tips with height comments", "[Tree]") {
    namu::TreeManip::SharedPtrVector trees;
    // Keeping original branch lengths here to make sure height comments are
    // used to override these
    std::string newick = "(1[&height=0]:0.3,2[&height=0.05]:0.3,(3[&height=0]:0.2,(4[&height=0.02]:0.1,5[&height=0]:0.1)[&height=0.1,height_index=0]:0.1)[&height=0.2,height_index=1]:0.1)[&height=0.3,height_index=2];";
    namu::TreeIO tree_io = namu::TreeIO();
    trees = tree_io.parse_from_newick(newick);
    REQUIRE(trees.size() == 1);
    namu::TreeManip::SharedPtr tm = trees.at(0);
    std::string new_str = tm->make_newick(3);
    std::string expected_new_str = "(1:0.300,2:0.250,(3:0.200,(4:0.080,5:0.100):0.100):0.100)";
    REQUIRE(new_str == expected_new_str);

    namu::Tree::SharedPtr tree = tm->get_tree();
    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);

    // 1.2 - 0.07
    REQUIRE(tm->calc_tree_length() == Approx(1.13));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    tm->scale_all_internal_heights(2.0);
    // 2.4 - 0.07
    REQUIRE(tm->calc_tree_length() == Approx(2.33));
    REQUIRE(tm->count_edges() == 7);
    REQUIRE(tm->count_nodes() == 8);

    REQUIRE(tree->num_leaves() == 5);
    REQUIRE(tree->num_internals() == 3);
    REQUIRE(tree->num_nodes() == 8);
    REQUIRE(tree->max_num_nodes() == 9);
}
