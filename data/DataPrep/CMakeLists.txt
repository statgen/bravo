cmake_minimum_required(VERSION 3.4)
project(DataPrep VERSION 0.0.1)
include(GNUInstallDirs)
add_compile_options(-Wall -Wpedantic)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(src)
add_subdirectory(test)

enable_testing()
add_test(NAME ${PROJECT_NAME}_test COMMAND testAll)

