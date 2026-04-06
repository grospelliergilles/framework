if (NOT TARGET mpi_target)
  return()
endif()

add_library(arcane_alina_mpi_example INTERFACE)
target_link_libraries(arcane_alina_mpi_example INTERFACE
  arcane_alina mpi_target ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_compile_definitions(arcane_alina_mpi_example INTERFACE "AMGCL_BLOCK_SIZES=${AMGCL_BLOCK_SIZES_SEQ}")

if (TARGET eigen_target)
  target_link_libraries(arcane_alina_mpi_example INTERFACE eigen_target)
endif()

if (TARGET Metis::metis)
  target_link_libraries(arcane_alina_mpi_example INTERFACE Metis::metis)
endif()

function(add_mpi_example example sources nb_core)
  add_executable(${example} ${sources})
  target_link_libraries(${example} PUBLIC arcane_alina_mpi_example)
  set(TEST_NAME "alina_sample_mpi${nb_core}_${example}")
  add_test(${TEST_NAME} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nb_core} ./${example})
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
add_mpi_example(mpi_complex mpi_complex.cc 2)

# Il faut un argument 'matrix' pour le test
add_mpi_example_no_test(schur_pc_mpi schur_pc_mpi.cc)

# Il faut un fichir de partition pour le test
add_mpi_example_no_test(solve_mm_mpi solve_mm_mpi.cc)

# Exception dans boost lors du test
add_mpi_example_no_test(cpr_mpi cpr_mpi.cc)

if (TARGET Metis::metis)
  add_executable(partition partition.cc)
  target_link_libraries(partition
    arcane_alina
    Metis::metis
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )
endif()

#target_link_libraries(call_mpi_lib PUBLIC arcane_alina_mpi_clib)

if (TARGET cuda_target)
  foreach(example runtime_sdd runtime_sdd_3d schur_pc_mpi mpi_solver)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/${example}.cu
      COPYONLY
    )
    cuda_add_executable(${example}_cuda
      ${CMAKE_CURRENT_BINARY_DIR}/${example}.cu
      OPTIONS -DSOLVER_BACKEND_CUDA)
    target_include_directories(${example}_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${example}_cuda PUBLIC arcane_alina_example cuda_target mpi_target)
  endforeach()
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_result
  ${CMAKE_CURRENT_BINARY_DIR}/plot_result
  COPYONLY
)
