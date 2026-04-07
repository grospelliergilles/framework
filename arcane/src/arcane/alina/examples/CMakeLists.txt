
set(ARCANE_ALINA_BLOCK_SIZES "3;4" CACHE STRING "Block sizes for block solvers in examples")
set(ARCANE_ALINA_BLOCK_SIZES_SEQ "")
foreach(B ${ARCANE_ALINA_BLOCK_SIZES})
  set(ARCANE_ALINA_BLOCK_SIZES_SEQ "${ARCANE_ALINA_BLOCK_SIZES_SEQ}(${B})")
endforeach()

function(add_arcane_alina_example_no_test EXAMPLE SOURCES)
  add_executable(${EXAMPLE} ${SOURCES})
  target_link_libraries(${EXAMPLE} PUBLIC arcane_alina Boost::program_options)
  target_include_directories(${EXAMPLE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../tests")
  target_compile_definitions(${EXAMPLE} PUBLIC "ARCANE_ALINA_BLOCK_SIZES=${ARCANE_ALINA_BLOCK_SIZES_SEQ}")
endfunction()

# NOTE:
# Certains tests utilisent OpenMP
# On fixe en dur le nombre de threads sinon on utilise tous
# les coeurs disponibles.
function(add_arcane_alina_example EXAMPLE SOURCES)
  add_arcane_alina_example_no_test(${EXAMPLE} ${SOURCES})
  #target_link_libraries(${EXAMPLE} PUBLIC arcane_alina_example)
  set(TEST_NAME "alina_sample_seq_${EXAMPLE}")
  add_test(NAME ${TEST_NAME} COMMAND ${EXAMPLE} ${ARGN})
  set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
endfunction()

add_arcane_alina_example_no_test(mm2bin mm2bin.cc)
add_arcane_alina_example_no_test(bin2mm bin2mm.cc)

add_arcane_alina_example(solver solver.cc)
add_arcane_alina_example(solver_complex solver_complex.cc)
add_arcane_alina_example(crs_builder adapter_crs_builder.cc)
add_arcane_alina_example(custom_adapter custom_adapter.cc)
add_arcane_alina_example(mixed_precision mixed_precision.cc)
add_arcane_alina_example(ublas adapter_ublas.cc)

set(SAMPLE_MATRIX_ARGS "--matrix" "${CMAKE_CURRENT_SOURCE_DIR}/../tests/bcsstk02.mtx")

add_arcane_alina_example_no_test(block_crs backend_block_crs.cc)
add_arcane_alina_example(cpr cpr.cc ${SAMPLE_MATRIX_ARGS})
add_arcane_alina_example(cpr_drs cpr_drs.cc ${SAMPLE_MATRIX_ARGS})
add_arcane_alina_example(ns_search ns_search.cc ${SAMPLE_MATRIX_ARGS})
add_arcane_alina_example_no_test(schur_pressure_correction schur_pressure_correction.cc)
add_arcane_alina_example_no_test(schurpc_mixed schurpc_mixed.cc)
add_arcane_alina_example_no_test(deflated_solver deflated_solver.cc)

target_link_libraries(ublas PUBLIC ${Boost_SERIALIZATION_LIBRARY})

if (TARGET eigen_target)
  add_arcane_alina_example(solver_eigen solver.cc)
  target_link_libraries(solver_eigen PUBLIC eigen_target)
  target_compile_definitions(solver_eigen PRIVATE SOLVER_BACKEND_EIGEN)

  add_arcane_alina_example_no_test(schurpc_mixed_eigen schurpc_mixed.cc)
  target_link_libraries(schurpc_mixed_eigen PUBLIC eigen_target)
  target_compile_definitions(schurpc_mixed_eigen PRIVATE BLOCK_TYPE_EIGEN)
endif()

add_arcane_alina_example(call_lib call_lib.cc)

if (TARGET cuda_target)
  enable_language(CUDA)
  function(add_cuda_example example)
    add_executable(${example}_cuda ${example}.cc)
    set_source_files_properties(${example}.cc PROPERTIES LANGUAGE CUDA)
    target_compile_definitions(${example}_cuda PUBLIC SOLVER_BACKEND_CUDA)
    target_include_directories(${example}_cuda PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../tests" )
    target_link_libraries(${example}_cuda PUBLIC arcane_alina cuda_target arccore_cuda_compile_flags Boost::program_options)
    add_test(NAME alina_${example}_cuda COMMAND ${example}_cuda ${SAMPLE_MATRIX_ARGS})
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
