#----------------------------------------------------------------------------
# Tests for Alina
#----------------------------------------------------------------------------
add_library(arcane_alina_test INTERFACE)

if (NOT Boost_USE_STATIC_LIBS)
  target_compile_definitions(arcane_alina_test INTERFACE BOOST_TEST_DYN_LINK)
endif ()

# ----------------------------------------------------------------------------

set(SOURCE_FILES
  test_qr.cc
  test_skyline_lu.cc
  test_complex_erf.cc
  test_solver_builtin.cc
  test_solver_complex.cc
  test_solver_block_crs.cc
  test_solver_ns_builtin.cc
  test_io.cc
  test_static_matrix.cc
)
if (TARGET eigen_target)
  list(APPEND SOURCE_FILES
    test_eigen_solver.cc
    test_solver_eigen.cc
    test_solver_ns_eigen.cc
  )
endif()

set(ACCELERATOR_FILES
)

arcane_add_component_test_executable(alina
  FILES ${SOURCE_FILES} ${ACCELERATOR_FILES}
)
target_include_directories(arcane_alina.tests INTERFACE ${Boost_INCLUDE_DIRS})

target_link_libraries(arcane_alina.tests PUBLIC
  arcane_alina
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  GTest::GTest
  GTest::Main
)

if (TARGET eigen_target)
  target_link_libraries(arcane_alina.tests PUBLIC eigen_target)
endif()


# NOTE:
# Certains tests utilisent OpenMP
# On fixe en dur le nombre de threads sinon on utilise tous
# les coeurs disponibles.

gtest_discover_tests(arcane_alina.tests
  DISCOVERY_TIMEOUT 30
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2"
)

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
