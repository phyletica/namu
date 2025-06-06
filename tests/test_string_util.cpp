#include <catch2/catch.hpp>
#include "namu/string_util.hpp"

SCENARIO("split provides Python-like splitting of strings", "[string_util]") {
    std::vector<std::string> elements;

    REQUIRE(elements.size() == 0);

    SECTION("splitting string into pre-allocated vector") {
        namu::str_util::split("Test_string__TO_Split", '_', elements);
        REQUIRE(elements.size() == 5);
        REQUIRE(elements[0] == "Test");
        REQUIRE(elements[1] == "string");
        REQUIRE(elements[2] == "");
        REQUIRE(elements[3] == "TO");
        REQUIRE(elements[4] == "Split");
    }

    SECTION("splitting string and returning elements in a new vector") {
        std::vector<std::string> words = namu::str_util::split("Test_string__TO_Split", '_');
        REQUIRE(words.size() == 5);
        REQUIRE(words[0] == "Test");
        REQUIRE(words[1] == "string");
        REQUIRE(words[2] == "");
        REQUIRE(words[3] == "TO");
        REQUIRE(words[4] == "Split");
    }

    SECTION("splitting string by tabs") {
        namu::str_util::split("Test\tstring\t\tTO\tSplit", '\t', elements);
        REQUIRE(elements.size() == 5);
        REQUIRE(elements[0] == "Test");
        REQUIRE(elements[1] == "string");
        REQUIRE(elements[2] == "");
        REQUIRE(elements[3] == "TO");
        REQUIRE(elements[4] == "Split");
    }
}

TEST_CASE("Testing pad_int", "[string_util]") {
    SECTION("Testing pad_int") {
        REQUIRE(namu::str_util::pad_int(0, 1) == "0");
        REQUIRE(namu::str_util::pad_int(0, 2) == "00");
        REQUIRE(namu::str_util::pad_int(0, 3) == "000");
        REQUIRE(namu::str_util::pad_int(1, 3) == "001");
        REQUIRE(namu::str_util::pad_int(11, 3) == "011");
        REQUIRE(namu::str_util::pad_int(111, 3) == "111");
        REQUIRE(namu::str_util::pad_int(1111, 3) == "1111");
        REQUIRE(namu::str_util::pad_int(1111, 4) == "1111");
        REQUIRE(namu::str_util::pad_int(1111, 5) == "01111");
    }
}

TEST_CASE("Testing get_indent", "[string_util]") {
    SECTION("Testing get_indent") {
        REQUIRE(namu::str_util::get_indent() == "    ");
        REQUIRE(namu::str_util::get_indent(0) == "");
        REQUIRE(namu::str_util::get_indent(1) == "    ");
        REQUIRE(namu::str_util::get_indent(2) == "        ");
        REQUIRE(namu::str_util::get_indent(3) == "            ");
    }
}

TEST_CASE("Testing lstrip", "[string_util]") {
    SECTION("Testing all white space") {
        std::string s = "\n\t   \t\n";
        std::string t = namu::str_util::lstrip(s);
        REQUIRE(s == "\n\t   \t\n");
        REQUIRE(t == "");
    }

    SECTION("Testing internal white space") {
        std::string s = "hello \t\n world";
        std::string t = namu::str_util::lstrip(s);
        REQUIRE(s == t);
    }
    
    SECTION("Testing left white space") {
        std::string s = " \t \nhello";
        std::string t = namu::str_util::lstrip(s);
        REQUIRE(s == " \t \nhello");
        REQUIRE(t == "hello");
    }

    SECTION("Testing right white space") {
        std::string s = "hello \t \n";
        std::string t = namu::str_util::lstrip(s);
        REQUIRE(s == t);
    }

    SECTION("Testing empty string") {
        std::string s = "";
        std::string t = namu::str_util::lstrip(s);
        REQUIRE(s == t);
    }
}

TEST_CASE("Testing rstrip", "[string_util]") {
    SECTION("Testing all white space") {
        std::string s = "\n\t   \t\n";
        std::string t = namu::str_util::rstrip(s);
        REQUIRE(s == "\n\t   \t\n");
        REQUIRE(t == "");
    }

    SECTION("Testing internal white space") {
        std::string s = "hello \t\n world";
        std::string t = namu::str_util::rstrip(s);
        REQUIRE(s == t);
    }
    
    SECTION("Testing left white space") {
        std::string s = " \t \nhello";
        std::string t = namu::str_util::rstrip(s);
        REQUIRE(s == t);
    }

    SECTION("Testing right white space") {
        std::string s = "hello \t \n";
        std::string t = namu::str_util::rstrip(s);
        REQUIRE(s == "hello \t \n");
        REQUIRE(t == "hello");
    }

    SECTION("Testing empty string") {
        std::string s = "";
        std::string t = namu::str_util::rstrip(s);
        REQUIRE(s == t);
    }
}

TEST_CASE("Testing strip", "[string_util]") {
    SECTION("Testing all white space") {
        std::string s = "\n\t   \t\n";
        std::string t = namu::str_util::strip(s);
        REQUIRE(s == "\n\t   \t\n");
        REQUIRE(t == "");
    }

    SECTION("Testing internal white space") {
        std::string s = "hello \t\n world";
        std::string t = namu::str_util::strip(s);
        REQUIRE(s == t);
    }

    SECTION("Testing left and right white space") {
        std::string s = "\n  \t \t hello \t \n";
        std::string t = namu::str_util::strip(s);
        REQUIRE(s == "\n  \t \t hello \t \n");
        REQUIRE(t == "hello");
    }

    SECTION("Testing empty string") {
        std::string s = "";
        std::string t = namu::str_util::strip(s);
        REQUIRE(s == t);
    }
}

TEST_CASE("Testing join", "[string_util]") {
    std::vector<std::string> strings {"1", "2", "3"};
    std::string s = namu::str_util::join(strings, "---");
    REQUIRE(s == "1---2---3");
}

TEST_CASE("Testing join on empty vector", "[string_util]") {
    std::vector<std::string> strings;
    std::string s = namu::str_util::join(strings, "---");
    REQUIRE(s == "");
}

TEST_CASE("Testing join on singleton vector", "[string_util]") {
    std::vector<std::string> strings {"1"};
    std::string s = namu::str_util::join(strings, "---");
    REQUIRE(s == "1");
}

TEST_CASE("Testing startswith", "[string_util]") {
    SECTION("Testing startswith") {
        REQUIRE(namu::str_util::startswith("amalie", "amalie") == true);
        REQUIRE(namu::str_util::startswith("amalie", "amali") == true);
        REQUIRE(namu::str_util::startswith("amalie", "a") == true);
        REQUIRE(namu::str_util::startswith("amaliee", "amalie") == true);
        REQUIRE(namu::str_util::startswith("amalie", "amaliee") == false);
        REQUIRE(namu::str_util::startswith("amalie", "A") == false);
        REQUIRE(namu::str_util::startswith("amalie", "eden") == false);
    }
}

TEST_CASE("Testing parse_map", "[string_util]") {
    SECTION("Testing parse_map") {
        std::string s = "a=1,b=2.00,foo=bar,key with spaces=value with spaces";
        std::map<std::string, std::string> expected_map;
        expected_map["a"] = "1";
        expected_map["b"] = "2.00";
        expected_map["foo"] = "bar";
        expected_map["key with spaces"] = "value with spaces";
        std::map<std::string, std::string> returned_map = namu::str_util::parse_map(s, ',', '=');
        REQUIRE(returned_map == expected_map);
    }
    SECTION("Testing parse_map with spaces") {
        std::string s = "  a = 1,   b     =    2.00  ,  foo = bar,key with spaces \t = \t  value with spaces    ";
        std::map<std::string, std::string> expected_map;
        expected_map["a"] = "1";
        expected_map["b"] = "2.00";
        expected_map["foo"] = "bar";
        expected_map["key with spaces"] = "value with spaces";
        std::map<std::string, std::string> returned_map = namu::str_util::parse_map(s, ',', '=');
        REQUIRE(returned_map == expected_map);
    }
    SECTION("Testing parse_map with adding to map") {
        std::vector<std::string> v = {
                "  a = 1   , b     =    2.00  ",
                " foo = bar",
                "\tkey with spaces   \t =  \t value with spaces  \t  "
        };
        std::map<std::string, std::string> expected_map;
        expected_map["a"] = "1";
        expected_map["b"] = "2.00";
        expected_map["foo"] = "bar";
        expected_map["key with spaces"] = "value with spaces";
        std::map<std::string, std::string> returned_map;
        for (auto s : v) {
            namu::str_util::parse_map(s, returned_map, ',', '=');
        }
        REQUIRE(returned_map == expected_map);
    }
}
