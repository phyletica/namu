#include <catch2/catch.hpp>
#include "phycoeval/string_util.hpp"
#include "phycoeval/error.hpp"

TEST_CASE("derived error classes can be thrown", "[error]") {

    std::string message = "Testing Error";

    SECTION("throwing PhycoevalError") {
        REQUIRE_THROWS_AS(throw phycoeval::PhycoevalError(message), phycoeval::PhycoevalError);
    }
}
