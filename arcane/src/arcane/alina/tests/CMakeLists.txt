#----------------------------------------------------------------------------
# The tests
#----------------------------------------------------------------------------
add_library(amgcl_test INTERFACE)

if (NOT Boost_USE_STATIC_LIBS)
  target_compile_definitions(amgcl_test INTERFACE BOOST_TEST_DYN_LINK)
endif ()

target_include_directories(amgcl_test INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(amgcl_test INTERFACE
  amgcl
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

# NOTE:
# Certains tests utilisent OpenMP
# On fixe en dur le nombre de threads sinon on utilise tous
# les coeurs disponibles.

function(add_amgcl_test TEST_EXEC_NAME TEST_SOURCE)
  add_executable(${TEST_EXEC_NAME} ${TEST_SOURCE})
  target_link_libraries(${TEST_EXEC_NAME} amgcl_test)
  set(TEST_NAME "alina_${TEST_EXEC_NAME}")
  add_test(${TEST_NAME} ${TEST_EXEC_NAME})
  set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
endfunction()

add_amgcl_test(test_skyline_lu test_skyline_lu.cc)
add_amgcl_test(test_complex_erf test_complex_erf.cc)
add_amgcl_test(test_qr test_qr.cc)
add_amgcl_test(test_solver_builtin test_solver_builtin.cc)
add_amgcl_test(test_solver_complex test_solver_complex.cc)
add_amgcl_test(test_solver_block_crs test_solver_block_crs.cc)
add_amgcl_test(test_solver_ns_builtin test_solver_ns_builtin.cc)
add_amgcl_test(test_io test_io.cc)

add_amgcl_test(test_static_matrix test_static_matrix.cc)
# target_compile_options(test_static_matrix PRIVATE
#   $<$<CXX_COMPILER_ID:GNU>:-std=c++0x>
#   $<$<CXX_COMPILER_ID:Clang>:-std=c++0x>
# )

if (TARGET eigen_target)
  add_amgcl_test(test_eigen_solver test_eigen_solver.cc)
  add_amgcl_test(test_solver_eigen test_solver_eigen.cc)
  add_amgcl_test(test_solver_ns_eigen test_solver_ns_eigen.cc)

  target_link_libraries(test_eigen_solver eigen_target)
  target_link_libraries(test_solver_eigen eigen_target)
  target_link_libraries(test_solver_ns_eigen eigen_target)
endif()

if (AMGCL_HAVE_PYTHON AND NOT WIN32)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_pyamgcl.py
    ${CMAKE_CURRENT_BINARY_DIR}/test_pyamgcl.py
    COPYONLY
  )
  add_test(test_pyamgcl test_pyamgcl.py)
  set_tests_properties(test_pyamgcl
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/..")
endif()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
