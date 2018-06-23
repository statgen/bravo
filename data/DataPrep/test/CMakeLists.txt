find_package(Threads REQUIRED)
find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(Z_LIB NAMES libz.a z HINT ${CGET_PREFIX}/lib)
find_library(BOOST_PROGRAM_OPTIONS_LIB boost_program_options HINT ${CGET_PREFIX}/lib)
find_library(CURL_LIB curl) # must be installed in system libraries
find_library(LZMA_LIB lzma) # must be installed in system libraries
find_library(LBZ2_LIB bz2) # must be installed in system libraries
find_library(GTEST_LIB gtest HINT ${CGET_PREFIX}/lib)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
    include_directories(${CGET_PREFIX}/include/htslib)
endif()

add_executable(testAll HistogramTest.cpp PercentilesTest.cpp Main_TestAll.cpp)
target_link_libraries(testAll Histogram Percentiles ${HTS_LIB} ${CURL_LIB} ${LZMA_LIB} ${LBZ2_LIB} ${Z_LIB} ${GTEST_LIB})
install(TARGETS testAll COMPONENT cli RUNTIME DESTINATION bin)
