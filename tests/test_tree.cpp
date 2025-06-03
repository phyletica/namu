#include <catch2/catch.hpp>
#include "phycoeval/tree.hpp"

TEST_CASE("Testing Tree bare constructor", "[Tree]") {
    phycoeval::Tree t = phycoeval::Tree();
    REQUIRE(t.num_nodes() == 0);
}
