# Unit tests. Each is a small Fortran program linked against the
# disorder physics code (disorder_core) and the no-op analysis; it
# exits non-zero if any of its checks fail.
add_library(disorder_test_support OBJECT
  test_utils.f90 ${PROJECT_SOURCE_DIR}/analysis/simple_analysis.f)
target_include_directories(disorder_test_support PUBLIC ${disorder_INCLUDES})
add_dependencies(disorder_test_support disorder_core)

function(disorder_add_test_program name)
  add_executable(${name} ${name}.f90)
  target_include_directories(${name} PUBLIC ${disorder_INCLUDES})
  target_link_libraries(${name} disorder_test_support disorder_core ${disorder_LIBRARIES})
endfunction()

foreach(t test_boosts test_phase_space test_disent_interface)
  disorder_add_test_program(${t})
  add_test(NAME unit_${t} COMMAND ${t})
  set_tests_properties(unit_${t} PROPERTIES LABELS "unit")
endforeach()

# The matrix element against LO reference expressions, once per process.
# The kinematic ranges only set up the HOPPET grids; the (x,Q) points
# are chosen inside the test.
disorder_add_test_program(test_matrix_element)
set(me_common -lo -toyQ0 2.0 -yorder 5 -lnlnQorder 4 -xmin 1e-3 -xmax 0.9 -Qmin 2.0)
set(me_cases
  "nc_photon_eminus|-NC"
  "nc_photon_eplus|-NC -positron"
  "nc_gammaZ_eminus|-NC -includeZ"
  "nc_gammaZ_eplus|-NC -includeZ -positron"
  "nc_Zonly_eminus|-NC -includeZ -Zonly"
  "nc_intonly_eplus|-NC -includeZ -intonly -positron"
  "nc_neutrino|-NC -neutrino"
  "nc_antineutrino|-NC -neutrino -positron"
  "cc_eminus|-noNC -CC"
  "cc_eplus|-noNC -CC -positron"
  "cc_neutrino|-noNC -CC -neutrino"
  "cc_antineutrino|-noNC -CC -neutrino -positron"
  "nc_cc_eminus|-NC -includeZ -CC")
# ... and beyond LO the combination of structure functions (photon only)
set(me_nlo_cases
  "nc_photon_nlo|-NC"
  "nc_photon_nnlo_eplus|-NC -positron")
foreach(case ${me_cases})
  string(REPLACE "|" ";" case ${case})
  list(GET case 0 name)
  list(GET case 1 args)
  separate_arguments(args)
  add_test(NAME unit_matrix_element_${name} COMMAND test_matrix_element ${me_common} ${args})
  set_tests_properties(unit_matrix_element_${name} PROPERTIES LABELS "unit")
endforeach()
string(REPLACE "-lo;" "" me_common_nlo "${me_common}")
foreach(case ${me_nlo_cases})
  string(REPLACE "|" ";" case ${case})
  list(GET case 0 name)
  list(GET case 1 args)
  separate_arguments(args)
  if(name MATCHES "nnlo")
    set(order -nnlo)
  else()
    set(order -nlo)
  endif()
  add_test(NAME unit_matrix_element_${name} COMMAND test_matrix_element ${order} ${me_common_nlo} ${args})
  set_tests_properties(unit_matrix_element_${name} PROPERTIES LABELS "unit")
endforeach()

# DISENT's three-parton matrix element for each process
disorder_add_test_program(test_matthr)
foreach(case ${me_cases})
  string(REPLACE "|" ";" case ${case})
  list(GET case 0 name)
  list(GET case 1 args)
  separate_arguments(args)
  add_test(NAME unit_matthr_${name} COMMAND test_matthr ${me_common} ${args})
  set_tests_properties(unit_matthr_${name} PROPERTIES LABELS "unit")
endforeach()

# Combinations that disorder must refuse to run. The guards in
# set_parameters use `stop`, which exits with status 0, so the tests
# look for the message instead.
set(guard_cases
  "p2b_cc|Cannot do CC in p2b beyond NLO yet|-p2b -CC -nnlo -toyQ0 2.0 -Q 10 -x 0.01"
  "p2b_n3lo|Cannot run p2b at N3LO yet|-p2b -n3lo -toyQ0 2.0 -Q 10 -x 0.01"
  "p2b_Z|Cannot do Z in p2b beyond NLO yet|-p2b -includeZ -nnlo -toyQ0 2.0 -Q 10 -x 0.01"
  "p2b_nu_nc|Cannot do Z in p2b beyond NLO yet|-p2b -NC -neutrino -nnlo -toyQ0 2.0 -Q 10 -x 0.01"
  "p2b_vnf|Cannot run p2b with variable flavour|-p2b -vnf -nnlo -toyQ0 2.0 -Q 10 -x 0.01"
  "p2b_pdfuncert|Cannot do pdf uncertainties with p2b|-p2b -nnlo -toyQ0 2.0 -Q 10 -x 0.01 -pdfuncert"
  "alphasuncert_alone|Need -pdfuncert to run with -alphasuncert|-nnlo -toyQ0 2.0 -Q 10 -x 0.01 -alphasuncert"
  "three_constraints|Only two variables can be constrained|-nnlo -toyQ0 2.0 -Q 10 -x 0.01 -y 0.1"
  "no_phase_space|No phase space available|-nnlo -toyQ0 2.0 -Q 400 -x 0.01")
foreach(case ${guard_cases})
  string(REPLACE "|" ";" case ${case})
  list(GET case 0 name)
  list(GET case 1 message)
  list(GET case 2 args)
  separate_arguments(args)
  add_test(NAME guard_${name} COMMAND disorder ${args})
  set_tests_properties(guard_${name} PROPERTIES LABELS "unit"
    PASS_REGULAR_EXPRESSION "${message}")
endforeach()
