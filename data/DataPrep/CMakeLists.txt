cmake_minimum_required(VERSION 3.1)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-O3")

project(SitesPreparation)

add_subdirectory(src)
add_subdirectory(test)

enable_testing()
add_test(NAME ${PROJECT_NAME}_test COMMAND testAll)

