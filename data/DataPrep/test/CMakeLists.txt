find_package(Threads REQUIRED)
find_library(Z_LIB NAMES libz.a z HINT ${CGET_PREFIX}/lib)
find_library(BOOST_PROGRAM_OPTIONS_LIB boost_program_options HINT ${CGET_PREFIX}/lib)
find_library(GTEST_LIB gtest HINT ${CGET_PREFIX}/lib)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

add_executable(testAll HistogramTest.cpp PercentilesTest.cpp Main_TestAll.cpp)
target_link_libraries(testAll Histogram Percentiles ${Z_LIB} ${GTEST_LIB})
install(TARGETS testAll COMPONENT cli RUNTIME DESTINATION bin)
