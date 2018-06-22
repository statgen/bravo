set(AUX_SRCS TypeSwitcher.cpp Histogram.cpp Percentiles.cpp Metric.cpp GzipWriter.cpp aux.cpp)

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
add_executable(CountAlleles CountAlleles.cpp ${AUX_SRCS})
add_executable(GTHistograms GTHistograms.cpp ${AUX_SRCS})
add_executable(INFOHistograms INFOHistograms.cpp ${AUX_SRCS})
add_executable(INFOPercentiles INFOPercentiles.cpp ${AUX_SRCS})
add_executable(RandomHetHom RandomHetHom.cpp ${AUX_SRCS})

target_link_libraries(CountAlleles ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(GTHistograms ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(INFOHistograms ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(INFOPercentiles ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})
target_link_libraries(RandomHetHom ${Z_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB} ${LBZ2_LIB} ${CURL_LIB} ${LZMA_LIB})

install(TARGETS CountAlleles GTHistograms INFOHistograms INFOPercentiles RandomHetHom COMPONENT cli RUNTIME DESTINATION bin)
