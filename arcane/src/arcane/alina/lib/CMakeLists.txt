add_library(arcane_alina_clib amgcl.cc)
target_link_libraries(arcane_alina_clib PUBLIC arcane_alina)
#set_target_properties(arcane_alina_clib PROPERTIES PREFIX "")
target_include_directories(arcane_alina_clib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (TARGET mpi_target)
  add_library(arcane_alina_mpi_clib amgcl_mpi.cc)
  #set_target_properties(arcane_alina_mpi_clib PROPERTIES PREFIX "")
  target_link_libraries(arcane_alina_mpi_clib PUBLIC arcane_alina_clib mpi_target)
endif()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
