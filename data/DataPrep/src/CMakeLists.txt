set(AUX_SRCS TypeSwitcher.cpp Histogram.cpp Percentiles.cpp aux.cpp)

find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(BOOST_PROGRAM_OPTIONS_LIB boost_program_options HINT ${CGET_PREFIX}/lib)
find_library(CURL_LIB curl) # must be installed in system libraries
find_library(CRYPTO_LIB crypto) # must be installed in system libraries 
find_library(LZMA_LIB lzma) # must be installed in system libraries
find_library(LBZ2_LIB bz2) # must be installed in system libraries
find_library(Z_LIB NAMES libz z) # must be installed in system libraries

find_package(Threads)

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

target_link_libraries(ComputeAlleleCounts ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${Z_LIB} ${CURL_LIB} ${LZMA_LIB} ${CRYPTO_LIB} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ComputeHistograms ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${Z_LIB}  ${CURL_LIB} ${LZMA_LIB} ${CRYPTO_LIB} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ComputeAlleleCountsAndHistograms ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${Z_LIB}  ${CURL_LIB} ${LZMA_LIB} ${CRYPTO_LIB} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ComputePercentiles ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${Z_LIB}  ${CURL_LIB} ${LZMA_LIB} ${CRYPTO_LIB} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(RandomHetHom ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${Z_LIB}  ${CURL_LIB} ${LZMA_LIB} ${CRYPTO_LIB} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS ComputeAlleleCounts ComputeHistograms ComputeAlleleCountsAndHistograms ComputePercentiles RandomHetHom COMPONENT cli RUNTIME DESTINATION bin)
