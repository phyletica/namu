#include <catch2/catch.hpp>
#include "namu/node.hpp"

TEST_CASE("Testing Node bare constructor", "[Node]") {
    namu::Node n = namu::Node();
    REQUIRE(n.is_root() == true);
    REQUIRE(n.is_leaf() == true);
    REQUIRE(n.has_parent() == false);
    REQUIRE(n.has_child() == false);
    REQUIRE(n.get_height() == 0.0);
    REQUIRE(n.get_edge_length() == 0.0);

    SECTION("Testing set_height") {
        REQUIRE(n.get_height() == 0.0);
        REQUIRE(n.get_edge_length() == 0.0);
        n.set_height(1.0);
        REQUIRE(n.get_height() == 1.0);
        REQUIRE(n.get_edge_length() == 0.0);
    }
}
