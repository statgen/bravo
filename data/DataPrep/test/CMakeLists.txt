find_package(Threads REQUIRED)
find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(BOOST_SYSTEM_LIB boost_system HINT ${CGET_PREFIX}/lib)
find_library(CURL_LIB curl) # must be installed in system libraries
find_library(CRYPTO_LIB crypto) # must be installed in system libraries
find_library(LZMA_LIB lzma) # must be installed in system libraries
find_library(LBZ2_LIB bz2) # must be installed in system libraries
find_library(Z_LIB NAMES libz z) # must be installed in system libraries
find_library(GTEST_LIB gtest HINT ${CGET_PREFIX}/lib)

find_package(Threads)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
    include_directories(${CGET_PREFIX}/include/htslib)
endif()

add_executable(testAll HistogramTest.cpp PercentilesTest.cpp ExecutableTest.cpp Main_TestAll.cpp)
target_link_libraries(testAll Histogram Percentiles ${BOOST_SYSTEM_LIB} ${HTS_LIB} ${CURL_LIB} ${LZMA_LIB} ${LBZ2_LIB} ${Z_LIB} ${GTEST_LIB} ${CRYPTO_LIB} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS testAll COMPONENT cli RUNTIME DESTINATION test)
install(FILES
            input.vcf.gz
            input.vcf.gz.tbi
            input.subset.vcf.gz
            input.subset.vcf.gz.tbi
            samples.txt
            allelecounts.output_ok.vcf.gz
            allelecounts.region.output_ok.vcf.gz
            histograms.output_ok.vcf.gz
            histograms.region.output_ok.vcf.gz
            allelecounts_and_histograms.output_ok.vcf.gz
        DESTINATION test)
