cmake_minimum_required(VERSION 2.6)

project(simplex-search)

add_definitions(
    -std=c++11
    )

enable_testing()

add_subdirectory(gtest-1.6.0)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_library(simplex_search
    simplex_search.cc
    )

add_executable(simplex_test simplex_search_test.cc)
target_link_libraries(simplex_test 
    simplex_search
    gtest
    gtest_main
    )
add_test(simplex_test simplex_test)
