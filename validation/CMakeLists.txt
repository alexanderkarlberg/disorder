# The validation matrix (configurations.txt), one ctest test per
# configuration, labelled "validation" plus its labels from the file
# ("fast" or "slow"). Each test runs disorder and compares its output
# with ref_runs/ (see run_validation.py).
#
# The reference histograms come from exclusive_lab_frame_analysis.f, so
# the tests are only registered for a build with that analysis.
if(NOT ANALYSIS STREQUAL "exclusive_lab_frame_analysis.f" OR NOT NEEDS_FASTJET)
  message(STATUS "Validation tests not registered: they need "
    "-DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f")
  return()
endif()
find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
  message(WARNING "Validation tests not registered: python3 not found")
  return()
endif()

# Where LHAPDF looks for PDF sets, to disable runs whose set is missing
execute_process(COMMAND ${LHAPDF_CONFIG} --datadir
  OUTPUT_VARIABLE lhapdf_datadir OUTPUT_STRIP_TRAILING_WHITESPACE)
set(lhapdf_paths ${lhapdf_datadir})
if(DEFINED ENV{LHAPDF_DATA_PATH})
  string(REPLACE ":" ";" env_paths "$ENV{LHAPDF_DATA_PATH}")
  list(APPEND lhapdf_paths ${env_paths})
endif()

file(STRINGS configurations.txt config_lines ENCODING UTF-8)
foreach(line ${config_lines})
  if(line MATCHES "^[ \t]*(#|$)")
    continue()
  endif()
  string(REPLACE "|" ";" fields "${line}")
  list(GET fields 0 labels)
  list(GET fields 1 prefix)
  list(GET fields 2 args)
  string(STRIP "${prefix}" prefix)
  separate_arguments(labels)
  string(REGEX REPLACE "_$" "" name "${prefix}")
  add_test(NAME validation_${name}
    COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/run_validation.py
            --disorder $<TARGET_FILE:disorder> --prefix ${prefix}
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/${name})
  set_tests_properties(validation_${name} PROPERTIES LABELS "validation;${labels}")

  if(args MATCHES "-pdf +([^ ]+)")
    set(pdfset ${CMAKE_MATCH_1})
    set(found FALSE)
    foreach(path ${lhapdf_paths})
      if(EXISTS ${path}/${pdfset})
        set(found TRUE)
      endif()
    endforeach()
    if(NOT found)
      message(WARNING "LHAPDF set ${pdfset} not found in ${lhapdf_paths}: "
        "validation_${name} is disabled")
      set_tests_properties(validation_${name} PROPERTIES DISABLED TRUE)
    endif()
  endif()
endforeach()
