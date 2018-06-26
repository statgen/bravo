set(AUX_SRCS TypeSwitcher.cpp Histogram.cpp Percentiles.cpp aux.cpp)

find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(Z_LIB NAMES libz.a z HINT ${CGET_PREFIX}/lib)
find_library(BOOST_PROGRAM_OPTIONS_LIB boost_program_options HINT ${CGET_PREFIX}/lib)
find_library(CURL_LIB curl) # must be installed in system libraries
find_library(LZMA_LIB lzma) # must be installed in system libraries
find_library(LBZ2_LIB bz2) # must be installed in system libraries

message("he ${CURL_LIB} ${LZMA_LIB} ${LBZ2_LIB}")

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
    include_directories(${CGET_PREFIX}/include/htslib)
endif()

add_library(Histogram Histogram.cpp)
add_library(Percentiles Percentiles.cpp ${AUX_SRCS})
add_executable(ComputeAlleleCounts ComputeAlleleCounts.cpp ${AUX_SRCS})
add_executable(ComputeHistograms ComputeHistograms.cpp ${AUX_SRCS})
add_executable(ComputeAlleleCountsAndHistograms ComputeAlleleCountsAndHistograms.cpp ${AUX_SRCS})
add_executable(ComputePercentiles ComputePercentiles.cpp ${AUX_SRCS})
add_executable(RandomHetHom RandomHetHom.cpp ${AUX_SRCS})

target_link_libraries(ComputeAlleleCounts ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(ComputeHistograms ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(ComputeAlleleCountsAndHistograms ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(ComputePercentiles ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(RandomHetHom ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})

install(TARGETS ComputeAlleleCounts ComputeHistograms ComputeAlleleCountsAndHistograms ComputePercentiles RandomHetHom COMPONENT cli RUNTIME DESTINATION bin)
