//
// Created by Joe Yesselman on 9/20/22.
//

#include "../common.hpp"
#include <iostream>
#include "../src/base/functions.h"

TEST_CASE("test ring mapper algorithm") {
    String filename("../src/unittest_resources/ring_mapper_results_reference.csv");
    std::ifstream ref_file(filename);
    //std::ifstream file_counter(filename);

    String line_ref;
    String line_test;

    SUBCASE("test original/reference mapper") {
        String filename_2("../src/unittest_resources/ring_mapper_results_py_algo.csv");
        std::ifstream test_file(filename_2);
        int counter = 0;
        Strings line_ref_vector;
        Strings line_test_vector;
        std::vector<double> line_ref_vec_double;
        line_ref_vec_double = {0, 0, 0};
        std::vector<double> line_test_vec_double;
        line_test_vec_double = {0, 0, 0};

        while (counter < 2478) {
            if (counter >= 1) {
                getline(ref_file, line_ref);
                line_ref_vector = split(line_ref, ",");

                getline(test_file, line_test);
                line_test_vector = split(line_test, ",");

                CHECK(line_ref_vec_double[0] == doctest::Approx(line_test_vec_double[0]));
                CHECK(line_ref_vec_double[1] == doctest::Approx(line_test_vec_double[1]));
                CHECK(line_ref_vec_double[2] == doctest::Approx(line_test_vec_double[2]));
            }
            counter++;
        }
    }
    SUBCASE("test basic mapper") {
        String filename_2("../src/unittest_resources/");
    }


    SUBCASE("test more advanced mapper") {
        String filename_2("insert test csv directory");
        std::ifstream test_file(filename_2);

    }
}