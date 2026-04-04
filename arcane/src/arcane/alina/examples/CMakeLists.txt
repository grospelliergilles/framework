add_library(amgcl_example INTERFACE)
target_include_directories(amgcl_example INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/../tests
  ${Boost_INCLUDE_DIRS}
)

target_link_libraries(amgcl_example INTERFACE
  amgcl
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

set(AMGCL_BLOCK_SIZES "3;4" CACHE STRING "Block sizes for block solvers in examples")
set(AMGCL_BLOCK_SIZES_SEQ "")
foreach(B ${AMGCL_BLOCK_SIZES})
  set(AMGCL_BLOCK_SIZES_SEQ "${AMGCL_BLOCK_SIZES_SEQ}(${B})")
endforeach()
target_compile_definitions(amgcl_example INTERFACE "AMGCL_BLOCK_SIZES=${AMGCL_BLOCK_SIZES_SEQ}")

function(add_amgcl_example_no_test EXAMPLE SOURCES)
  add_executable(${EXAMPLE} ${SOURCES})
  target_link_libraries(${EXAMPLE} PUBLIC amgcl_example)
endfunction()

# NOTE:
# Certains tests utilisent OpenMP
# On fixe en dur le nombre de threads sinon on utilise tous
# les coeurs disponibles.
function(add_amgcl_example EXAMPLE SOURCES)
  add_executable(${EXAMPLE} ${SOURCES})
  target_link_libraries(${EXAMPLE} PUBLIC amgcl_example)
  set(TEST_NAME "alina_sample_seq_${EXAMPLE}")
  add_test(${TEST_NAME} ${EXAMPLE})
  set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
endfunction()

add_amgcl_example_no_test(mm2bin mm2bin.cc)
add_amgcl_example_no_test(bin2mm bin2mm.cc)

add_amgcl_example(solver solver.cc)
add_amgcl_example(solver_complex solver_complex.cc)
add_amgcl_example(crs_builder adapter_crs_builder.cc)
add_amgcl_example(custom_adapter custom_adapter.cc)
add_amgcl_example(mixed_precision mixed_precision.cc)
add_amgcl_example(ublas adapter_ublas.cc)

add_amgcl_example_no_test(block_crs block_crs.cc)
add_amgcl_example_no_test(cpr cpr.cc)
add_amgcl_example_no_test(cpr_drs cpr_drs.cc)
add_amgcl_example_no_test(schur_pressure_correction schur_pressure_correction.cc)
add_amgcl_example_no_test(schurpc_mixed schurpc_mixed.cc)
add_amgcl_example_no_test(deflated_solver deflated_solver.cc)
add_amgcl_example_no_test(ns_search ns_search.cc)

target_link_libraries(ublas PUBLIC ${Boost_SERIALIZATION_LIBRARY})

if (TARGET eigen_target)
  add_amgcl_example(solver_eigen solver.cc)
  target_link_libraries(solver_eigen PUBLIC eigen_target)
  target_compile_definitions(solver_eigen PRIVATE SOLVER_BACKEND_EIGEN)

  add_amgcl_example_no_test(schurpc_mixed_eigen schurpc_mixed.cc)
  target_link_libraries(schurpc_mixed_eigen PUBLIC eigen_target)
  target_compile_definitions(schurpc_mixed_eigen PRIVATE BLOCK_TYPE_EIGEN)
endif()

add_amgcl_example(call_lib call_lib.cc)
target_link_libraries(call_lib PUBLIC libamgcl)


if (TARGET cuda_target)
  function(add_cuda_example example)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/${example}.cu
      COPYONLY
    )

    cuda_add_executable(${example}_cuda ${CMAKE_CURRENT_BINARY_DIR}/${example}.cu
      OPTIONS -DSOLVER_BACKEND_CUDA)
    target_link_libraries(${example}_cuda amgcl_example cuda_target)
  endfunction()

  add_cuda_example(solver)
  add_cuda_example(schur_pressure_correction)
  add_cuda_example(cpr_drs)
endif()

add_subdirectory(mpi)

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
