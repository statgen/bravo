set(AUX_SRCS TypeSwitcher.cpp Histogram.cpp Percentiles.cpp Metric.cpp GzipWriter.cpp aux.cpp)

find_path(BOOST_PROGRAM_OPTIONS_INCL NAMES program_options.hpp PATH_SUFFIXES include include/boost)
if(NOT BOOST_PROGRAM_OPTIONS_INCL)
   message(SEND_ERROR "Boost.Program_options header was not found!")
else(NOT BOOST_PROGRAM_OPTIONS_INCL)
   message(STATUS "Boost.Program_options header directory: " ${BOOST_PROGRAM_OPTIONS_INCL})
endif(NOT BOOST_PROGRAM_OPTIONS_INCL)

find_library(BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options PATH_SUFFIXES lib)
if(NOT BOOST_PROGRAM_OPTIONS_LIB)
   message(SEND_ERROR "Boost.Program_options library was not found!")
else(NOT BOOST_PROGRAM_OPTIONS_LIB)
   message(STATUS "Boost.Program_options library: " ${BOOST_PROGRAM_OPTIONS_LIB})
endif(NOT BOOST_PROGRAM_OPTIONS_LIB)

find_path(HTS_INCL NAMES hts.h vcf.h synced_bcf_reader.h PATH_SUFFIXES include include/htslib)
if(NOT HTS_INCL)
   message(SEND_ERROR "HTSlib headers were not found!")
else(NOT HTS_INCL)
   message(STATUS "HTSlib headers directory: " ${HTS_INCL})
endif(NOT HTS_INCL)

find_library(HTS_LIB NAMES hts PATH_SUFFIXES lib)
if(NOT HTS_LIB)
   message(SEND_ERROR "HTSlib library was not found!")
else(NOT HTS_LIB)
   message(STATUS "HTSlib library: " ${HTS_LIB})
endif(NOT HTS_LIB)

find_package(ZLIB REQUIRED)

message(STATUS "ZLIB version: " ${ZLIB_VERSION_STRING})
message(STATUS "ZLIB header directory: " ${ZLIB_INCLUDE_DIRS})
message(STATUS "ZLIB library: " ${ZLIB_LIBRARIES})

include_directories(${ZLIB_INCLUDE_DIRS} ${HTS_INCL} ${BOOST_PROGRAM_OPTIONS_INCL})

add_library(Histogram Histogram.cpp)
add_library(Percentiles Percentiles.cpp ${AUX_SRCS})
add_executable(CountAlleles CountAlleles.cpp ${AUX_SRCS})
add_executable(GTHistograms GTHistograms.cpp ${AUX_SRCS})
add_executable(INFOHistograms INFOHistograms.cpp ${AUX_SRCS})
add_executable(INFOPercentiles INFOPercentiles.cpp ${AUX_SRCS})
add_executable(RandomHetHom RandomHetHom.cpp ${AUX_SRCS})

target_link_libraries(CountAlleles ${ZLIB_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB})
target_link_libraries(GTHistograms ${ZLIB_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB})
target_link_libraries(INFOHistograms ${ZLIB_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB})
target_link_libraries(INFOPercentiles ${ZLIB_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB})
target_link_libraries(RandomHetHom ${ZLIB_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIB} ${HTS_LIB})

install(TARGETS CountAlleles GTHistograms INFOHistograms INFOPercentiles RandomHetHom DESTINATION ${CMAKE_BINARY_DIR}/bin)
