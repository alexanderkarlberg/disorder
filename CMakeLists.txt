cmake_minimum_required(VERSION 3.7.0)
include(CMakePrintHelpers)
project(disorder VERSION 1.0.0 LANGUAGES CXX Fortran)

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

set(CMAKE_Fortran_MODULE_DIRECTORY modules)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -ffixed-line-length-132 -fPIC")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra -fPIC")

set(CMAKE_BUILD_TYPE Release)
set(prefix ${CMAKE_INSTALL_PREFIX})
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (disorder_INCLUDES )
set (disorder_SOURCES )
set (disorder_ANALYSIS )
set (disorder_LIBRARIES )

# Main sources
file(GLOB_RECURSE disorder_SOURCES src/* analysis/pwhg_bookhist-multi.*)

# THE ANALYSIS TO BE COMPILED
file(GLOB_RECURSE disorder_ANALYSIS analysis/simple_analysis.f)
if(ANALYSIS)
  file(GLOB_RECURSE disorder_ANALYSIS analysis/${ANALYSIS})
endif()
cmake_print_variables(disorder_ANALYSIS)

list(APPEND disorder_SOURCES ${disorder_ANALYSIS})

##################################################################################
# Needs hoppet
find_file(HOPPET_CONFIG hoppet-config)
cmake_print_variables(HOPPET_CONFIG)
execute_process(
  COMMAND bash -c "${HOPPET_CONFIG} --prefix"	
  OUTPUT_VARIABLE HOPPET_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE ret
  )

if(NOT ret EQUAL "0")
    message( FATAL_ERROR "Could not find hoppet-config in ${HOPPET_CONFIG}")
endif()

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

##################################################################################
#Needs LHAPDF
find_file(LHAPDF_CONFIG lhapdf-config)
cmake_print_variables(LHAPDF_CONFIG)
execute_process(
  COMMAND bash -c "${LHAPDF_CONFIG} --incdir"	
  OUTPUT_VARIABLE LHAPDF_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE ret
  )

if(NOT ret EQUAL "0")
    message( FATAL_ERROR "Could not find hoppet-config in ${LHAPDF_CONFIG}")
endif()

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

##################################################################################
#Check if we need fastjet
option(NEEDS_FASTJET "To link fastjet libraries" OFF)
cmake_print_variables(NEEDS_FASTJET)
if (NEEDS_FASTJET)
  file(GLOB_RECURSE fastjet_SOURCE  analysis/fastjetfortran.cc)
  list(APPEND disorder_SOURCES ${fastjet_SOURCE})
  find_file(FASTJET_CONFIG fastjet-config)
  cmake_print_variables(FASTJET_CONFIG)
  execute_process(
    COMMAND bash -c "${FASTJET_CONFIG} --prefix"	
    OUTPUT_VARIABLE FASTJET_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE ret
    )
  
  if(NOT ret EQUAL "0")
    message( FATAL_ERROR "Could not find hoppet-config in ${FASTJET_CONFIG}")
  endif()
  
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
endif()

cmake_print_variables(disorder_INCLUDES)
cmake_print_variables(disorder_LIBRARIES)

add_executable(disorder src/disorder.f90)
add_executable(mergedata aux/mergedata.f)
add_executable(getpdfuncert aux/getpdfuncert.f)

target_sources(disorder PUBLIC ${disorder_SOURCES})
target_include_directories(disorder PUBLIC ${disorder_INCLUDES})
target_link_libraries(disorder ${disorder_LIBRARIES})
target_link_libraries(getpdfuncert ${LHAPDF_LIBRARY})

install(TARGETS disorder RUNTIME DESTINATION bin)


