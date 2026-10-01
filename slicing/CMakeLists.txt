# N-jettiness slicing studies (branch 2026-10-tau-slicing): SCET
# ingredients and the NLO DIS 2+1 test against DISENT.
add_library(slicing_scet OBJECT mod_slicing_scet.f90)
target_include_directories(slicing_scet PUBLIC ${disorder_INCLUDES})
add_dependencies(slicing_scet disorder_core)

foreach(prog test_scet tau2_nlo)
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${prog}.f90)
    add_executable(${prog} ${prog}.f90 ${PROJECT_SOURCE_DIR}/analysis/simple_analysis.f)
    target_include_directories(${prog} PUBLIC ${disorder_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(${prog} slicing_scet disorder_core ${disorder_LIBRARIES})
  endif()
endforeach()
