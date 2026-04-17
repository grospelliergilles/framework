
set(ARCCORE_ALINA_BLOCK_SIZES "3;4" CACHE STRING "Block sizes for block solvers in examples")
set(ARCCORE_ALINA_BLOCK_SIZES_SEQ "")
foreach(B ${ARCCORE_ALINA_BLOCK_SIZES})
  set(ARCCORE_ALINA_BLOCK_SIZES_SEQ "${ARCCORE_ALINA_BLOCK_SIZES_SEQ}(${B})")
endforeach()

# Add a library for common samples functions and types
add_library(arcane_alina_samples_lib
  AlinaSamplesCommon.h
  AlinaSamplesCommon.cc
)
arcane_target_set_standard_path(arcane_alina_samples_lib)
target_link_libraries(arcane_alina_samples_lib PRIVATE arcane_launcher arcane_full)
target_link_libraries(arcane_alina_samples_lib PUBLIC arccore_alina)
#target_compile_definitions(arcane_alina_samples_lib PRIVATE ARCANE_COMPONENT_arcane_alina)
target_include_directories(arcane_alina_samples_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

set(ALINA_SAMPLE_BUILD_DIR "${CMAKE_BINARY_DIR}/alina_samples")
function(add_arcane_alina_example_no_test EXAMPLE SOURCES)
  add_executable(${EXAMPLE} ${SOURCES})
  set_target_properties(${EXAMPLE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${ALINA_SAMPLE_BUILD_DIR}")
  target_link_libraries(${EXAMPLE} PUBLIC arcane_alina_samples_lib Boost::program_options)
  target_include_directories(${EXAMPLE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../tests")
  target_compile_definitions(${EXAMPLE} PUBLIC "ARCCORE_ALINA_BLOCK_SIZES=${ARCCORE_ALINA_BLOCK_SIZES_SEQ}")
endfunction()

# NOTE:
# Certains tests utilisent OpenMP
# On fixe en dur le nombre de threads sinon on utilise tous
# les coeurs disponibles par défaut.
function(add_arcane_alina_example EXAMPLE SOURCES)
  add_arcane_alina_example_no_test(${EXAMPLE} ${SOURCES})
  set(TEST_NAME "alina_sample_seq_${EXAMPLE}")
  add_test(NAME ${TEST_NAME}
    COMMAND ${EXAMPLE} ${ARGN}
    WORKING_DIRECTORY "${ALINA_SAMPLE_BUILD_DIR}"
  )
  set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
endfunction()

add_arcane_alina_example_no_test(mm2bin mm2bin.cc)
add_arcane_alina_example_no_test(bin2mm bin2mm.cc)

arcane_find_package(Hypre)

add_arcane_alina_example(basic_solver basic_solver.cc)
target_sources(basic_solver PUBLIC HypreComparer.cc)
arcane_add_arccon_packages(basic_solver PRIVATE Hypre)

add_arcane_alina_example(solver solver.cc)
target_sources(solver PUBLIC HypreComparer.cc)
arcane_add_arccon_packages(solver PRIVATE Hypre)

add_arcane_alina_example(solver_complex solver_complex.cc)
add_arcane_alina_example(crs_builder adapter_crs_builder.cc)
add_arcane_alina_example(custom_adapter custom_adapter.cc)
add_arcane_alina_example(mixed_precision mixed_precision.cc)

set(SAMPLE_MATRIX_ARGS "--matrix" "${CMAKE_CURRENT_SOURCE_DIR}/../tests/bcsstk02.mtx")

add_arcane_alina_example_no_test(block_crs backend_block_crs.cc)
add_arcane_alina_example(cpr cpr.cc ${SAMPLE_MATRIX_ARGS})
add_arcane_alina_example(cpr_drs cpr_drs.cc ${SAMPLE_MATRIX_ARGS})
add_arcane_alina_example(ns_search ns_search.cc ${SAMPLE_MATRIX_ARGS})
add_arcane_alina_example_no_test(schur_pressure_correction schur_pressure_correction.cc)
add_arcane_alina_example_no_test(schurpc_mixed schurpc_mixed.cc)
add_arcane_alina_example_no_test(deflated_solver DeflatedSolver.cc)
add_arcane_alina_example(poisson2d poisson2d.cc)

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
    target_link_libraries(${example}_cuda PUBLIC arccore_alina cuda_target arccore_cuda_compile_flags Boost::program_options)
    add_test(NAME alina_${example}_cuda COMMAND ${example}_cuda ${SAMPLE_MATRIX_ARGS})
  endfunction()

  add_cuda_example(solver)
  add_cuda_example(schur_pressure_correction)
  add_cuda_example(cpr_drs)
endif()

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Handle samples using MPI

add_library(arcane_alina_mpi_example INTERFACE)
target_link_libraries(arcane_alina_mpi_example INTERFACE
  arccore_alina ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_compile_definitions(arcane_alina_mpi_example INTERFACE "ARCCORE_ALINA_BLOCK_SIZES=${ARCCORE_ALINA_BLOCK_SIZES_SEQ}")

if (TARGET eigen_target)
  target_link_libraries(arcane_alina_mpi_example INTERFACE eigen_target)
endif()

function(add_mpi_example example sources nb_core)
  add_executable(${example} ${sources})
  set_target_properties(${example} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${ALINA_SAMPLE_BUILD_DIR}")
  target_link_libraries(${example} PUBLIC arcane_alina_samples_lib arcane_alina_mpi_example)
  set(TEST_NAME "alina_sample_mpi${nb_core}_${example}")
  add_test(NAME ${TEST_NAME}
    WORKING_DIRECTORY "${ALINA_SAMPLE_BUILD_DIR}"
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nb_core} ./${example}
  )
  set_property(TEST ${TEST_NAME} APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
endfunction()

function(add_mpi_example_no_test example sources)
  add_executable(${example} ${sources})
  target_link_libraries(${example} PUBLIC arcane_alina_mpi_example)
endfunction()

# Il faut mettre 1 seul coeur pour certains tests
# car ils utilisent openmp et sinon on se retrouve avec trop
# de processus

add_mpi_example(check_direct check_direct.cc 4)
add_mpi_example(runtime_sdd runtime_sdd.cc 2)
add_mpi_example(runtime_sdd_3d runtime_sdd_3d.cc 2)
add_mpi_example(runtime_bp runtime_bp.cc 2)
add_mpi_example(call_mpi_lib call_mpi_lib.cc 2)
add_mpi_example(test_transpose test_transpose.cc 4)
add_mpi_example(test_spmm test_spmm.cc 4)
add_mpi_example(spmm_scaling spmm_scaling.cc 4)
add_mpi_example(mpi_solver mpi_solver.cc 2)
add_mpi_example(mpi_basic_solver mpi_basic_solver.cc 2)
add_mpi_example(mpi_complex mpi_complex.cc 2)

# Il faut un argument 'matrix' pour le test
add_mpi_example_no_test(schur_pc_mpi schur_pc_mpi.cc)

# Il faut un fichir de partition pour le test
add_mpi_example_no_test(solve_mm_mpi solve_mm_mpi.cc)

# Exception dans boost lors du test
add_mpi_example_no_test(cpr_mpi cpr_mpi.cc)

#if(Parmetis_FOUND)
#  message(STATUS "Found Parmetis for Alina")
#  target_compile_definitions(arccore_alina PUBLIC ARCCORE_ALINA_HAVE_PARMETIS ARCCORE_ALINA_HAVE_PARMETIS)
#  target_link_libraries(arccore_alina PUBLIC arccon::Parmetis)
#endif()

if (TARGET arccon::Parmetis)
  add_mpi_example_no_test(partition partition.cc 2)
  target_link_libraries(partition PUBLIC arccon::Parmetis ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()

if (TARGET cuda_target)
  foreach(example runtime_sdd runtime_sdd_3d schur_pc_mpi mpi_solver)
    add_executable(${example}_cuda ${example}.cc)
    set_source_files_properties(${example}.cc PROPERTIES LANGUAGE CUDA)
    target_compile_definitions(${example}_cuda PUBLIC SOLVER_BACKEND_CUDA)
    target_link_libraries(${example}_cuda PUBLIC
      arcane_alina
      cuda_target
      arccore_cuda_compile_flags
      Boost::program_options
    )
  endforeach()
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_result
  ${CMAKE_CURRENT_BINARY_DIR}/plot_result
  COPYONLY
)

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
