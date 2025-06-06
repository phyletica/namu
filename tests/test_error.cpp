#include <catch2/catch.hpp>
#include "namu/string_util.hpp"
#include "namu/error.hpp"

TEST_CASE("derived error classes can be thrown", "[error]") {

    std::string message = "Testing Error";

    SECTION("throwing NamuX") {
        REQUIRE_THROWS_AS(throw namu::NamuX(message), namu::NamuX);
    }
}
