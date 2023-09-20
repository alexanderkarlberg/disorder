cmake_minimum_required(VERSION 3.10)
include(CMakePrintHelpers)
project(disorder VERSION 1.0.0 LANGUAGES CXX Fortran)

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -ffixed-line-length-132 -fPIC -mcmodel=large")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra -fPIC")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_Fortran_MODULE_DIRECTORY obj)
#SET(PROJECT_SOURCE_DIRECTORY src)

set (disorder_INCLUDES )
set (disorder_SOURCES )
set (disorder_ANALYSIS )
set (disorder_LIBRARIES )

# Main sources
file(GLOB_RECURSE  disorder_SOURCES src/* analysis/pwhg_bookhist-multi.*)
# THE ANALYSIS TO BE COMPILED
file(GLOB_RECURSE disorder_ANALYSIS analysis/cut_Ecur.f analysis/fastjetfortran.cc)
cmake_print_variables(disorder_ANALYSIS)

list(APPEND disorder_SOURCES ${disorder_ANALYSIS})

# Needs hoppet
find_file(HOPPET_CONFIG hoppet-config)
cmake_print_variables(HOPPET_CONFIG)
execute_process(
  COMMAND bash -c "${HOPPET_CONFIG} --prefix"	
  OUTPUT_VARIABLE HOPPET_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
set(HOPPET_INCLUDE_DIR "${HOPPET_INCLUDE_DIR}/include/hoppet")
cmake_print_variables(HOPPET_INCLUDE_DIR)

execute_process(
  COMMAND bash -c "${HOPPET_CONFIG} --libs"
  OUTPUT_VARIABLE HOPPET_LIBRARY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
cmake_print_variables(HOPPET_LIBRARY)

list(APPEND disorder_INCLUDES ${HOPPET_INCLUDE_DIR})
list(APPEND disorder_LIBRARIES ${HOPPET_LIBRARY})

#Needs LHAPDF
find_file(LHAPDF_CONFIG lhapdf-config)
cmake_print_variables(LHAPDF_CONFIG)
execute_process(
  COMMAND bash -c "${LHAPDF_CONFIG} --incdir"	
  OUTPUT_VARIABLE LHAPDF_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
set(LHAPDF_INCLUDE_DIR "${LHAPDF_INCLUDE_DIR}")
cmake_print_variables(LHAPDF_INCLUDE_DIR)

execute_process(
  COMMAND bash -c "${LHAPDF_CONFIG} --libs"
  OUTPUT_VARIABLE LHAPDF_LIBRARY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
cmake_print_variables(LHAPDF_LIBRARY)

list(APPEND disorder_INCLUDES ${LHAPDF_INCLUDE_DIR})
list(APPEND disorder_LIBRARIES ${LHAPDF_LIBRARY})

#Needs fastjet
find_file(FASTJET_CONFIG fastjet-config)
cmake_print_variables(FASTJET_CONFIG)
execute_process(
  COMMAND bash -c "${FASTJET_CONFIG} --prefix"	
  OUTPUT_VARIABLE FASTJET_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
set(FASTJET_INCLUDE_DIR "${FASTJET_INCLUDE_DIR}/include/fastjet")
cmake_print_variables(FASTJET_INCLUDE_DIR)

execute_process(
  COMMAND bash -c "${FASTJET_CONFIG} --libs --plugins"
  OUTPUT_VARIABLE FASTJET_LIBRARY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
cmake_print_variables(FASTJET_LIBRARY)

list(APPEND disorder_INCLUDES ${FASTJET_INCLUDE_DIR})
list(APPEND disorder_LIBRARIES ${FASTJET_LIBRARY})

cmake_print_variables(disorder_INCLUDES)
cmake_print_variables(disorder_LIBRARIES)

add_executable(disorder src/disorder.f90)
add_executable(mergedata aux/mergedata.f)
add_executable(getpdfuncert aux/getpdfuncert.f)

target_sources(disorder PUBLIC ${disorder_SOURCES})
target_include_directories(disorder PUBLIC ${disorder_INCLUDES})
target_link_libraries(disorder ${disorder_LIBRARIES})
target_link_libraries(getpdfuncert ${LHAPDF_LIBRARY})


