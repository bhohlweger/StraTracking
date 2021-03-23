cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(binFolder "${CMAKE_CURRENT_BINARY_DIR}/install")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${binFolder}/lib)
project(StraTrack CXX C)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(SRCS
	HarryPlotter.cc
)

set(HEADERS
	HarryPlotter.hh	
)

add_library(StraTrack_Helpers STATIC ${SRCS})

set_target_properties(StraTrack_Helpers PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(StraTrack_Helpers ${ROOT_LIBRARIES})

add_executable(ResolutionPVDCA  ${PROJECT_SOURCE_DIR}/ResolutionPVDCA.cxx)
target_link_libraries(ResolutionPVDCA StraTrack_Helpers)

add_executable(AddedHitPerformance ${PROJECT_SOURCE_DIR}/AddedHitPerformance.cxx)
target_link_libraries(AddedHitPerformance StraTrack_Helpers)

add_executable(ResolutionDV ${PROJECT_SOURCE_DIR}/ResolutionDV.cxx)
target_link_libraries(ResolutionDV StraTrack_Helpers)

add_executable(DCAProdCuts ${PROJECT_SOURCE_DIR}/DCAProdCuts.cxx)
target_link_libraries(DCAProdCuts StraTrack_Helpers)

add_executable(DCADaugCuts ${PROJECT_SOURCE_DIR}/DCADaugCuts.cxx)
target_link_libraries(DCADaugCuts StraTrack_Helpers)

add_executable(ResolutionDCADaug ${PROJECT_SOURCE_DIR}/ResolutionDCADaug.cxx)
target_link_libraries(ResolutionDCADaug StraTrack_Helpers)

add_executable(DCADaugProdCuts ${PROJECT_SOURCE_DIR}/DCADaugProdCuts.cxx)
target_link_libraries(DCADaugProdCuts StraTrack_Helpers)

