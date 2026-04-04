add_executable(poisson3Db poisson3Db.cc)
target_link_libraries(poisson3Db amgcl)

add_executable(poisson3Db_mixed poisson3Db.cc)
target_link_libraries(poisson3Db_mixed amgcl)
target_compile_definitions(poisson3Db_mixed PRIVATE MIXED_PRECISION)

if (TARGET cuda_target)
  cuda_add_executable(poisson3Db_cuda poisson3Db_cuda.cu)
  target_link_libraries(poisson3Db_cuda amgcl cuda_target)
endif()

if(TARGET mpi_target)
  add_executable(poisson3Db_mpi poisson3Db_mpi.cc)
  target_link_libraries(poisson3Db_mpi amgcl mpi_target)

  if (TARGET Metis::metis)
    target_link_libraries(poisson3Db_mpi Metis::metis)
  endif ()
endif()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
