include (srcs.cmake)

arcane_add_library(arcane_alina
  INPUT_PATH ${Arcane_SOURCE_DIR}/src
  RELATIVE_PATH arcane/alina
  FILES ${ARCANE_SOURCES}
)

#----------------------------------------------------------------------------
# Common compile options and definitions
#----------------------------------------------------------------------------
# target_compile_options(arcane_alina PUBLIC
#   # Compiler is GNU (g++):
#   $<$<CXX_COMPILER_ID:GNU>:$<BUILD_INTERFACE:-Wno-conversion;-Wno-float-conversion>>
#   # Compiler is Clang:
#   $<$<CXX_COMPILER_ID:Clang>:$<BUILD_INTERFACE:-Wno-conversion;-Wno-float-conversion>>
#   # Compiler is MSVC:
#   #$<$<CXX_COMPILER_ID:MSVC>:/bigobj>
#   #$<$<CXX_COMPILER_ID:MSVC>:/wd4715>
# )

target_compile_definitions(arcane_alina PRIVATE ARCANE_COMPONENT_arcane_alina)

target_link_libraries(arcane_alina PUBLIC Arcane::arcane_utils)
arcane_register_library(arcane_alina OPTIONAL)

if ((CMAKE_CXX_COMPILER_ID STREQUAL GNU) OR (CMAKE_CXX_COMPILER_ID STREQUAL Clang))
  target_compile_options(arcane_alina PUBLIC -Wno-conversion -Wno-float-conversion)
endif()

#----------------------------------------------------------------------------
# Find Boost
#----------------------------------------------------------------------------
option(Boost_USE_STATIC_LIBS "Use static versions of Boost libraries" OFF)
if (WIN32)
  set(Boost_USE_STATIC_LIBS ON)
endif ()

find_package(Boost REQUIRED COMPONENTS
  headers
  program_options
  serialization
  unit_test_framework
)

#----------------------------------------------------------------------------
# Builtin backend
#----------------------------------------------------------------------------
#add_library(amgcl INTERFACE)
#add_library(amgcl::amgcl ALIAS amgcl)
#target_sources(amgcl INTERFACE ${ARCANE_SOURCES})
#target_include_directories(amgcl INTERFACE ${Arcane_SOURCE_DIR}/src)
  
# target_compile_features(amgcl INTERFACE
#   cxx_auto_type
#   cxx_range_for
#   cxx_rvalue_references
#   cxx_right_angle_brackets
#   cxx_static_assert
# )

if (Boost_FOUND)
  message(STATUS "Found boost")
  target_include_directories(arcane_alina SYSTEM PUBLIC
    ${Boost_INCLUDE_DIRS}
  )
  target_link_libraries(arcane_alina PUBLIC Boost::headers)
else()
  target_compile_definitions(arcane_alina INTERFACE ARCANE_ALINA_NO_BOOST)
endif()

#target_include_directories(amgcl INTERFACE
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#  $<INSTALL_INTERFACE:include>
#)

find_package(OpenMP)
if (OPENMP_FOUND)
  target_compile_options(arcane_alina PUBLIC ${OpenMP_CXX_FLAGS})
  target_link_libraries(arcane_alina PUBLIC
    $<$<CXX_COMPILER_ID:GNU>:${OpenMP_CXX_FLAGS}>
    $<$<CXX_COMPILER_ID:Clang>:${OpenMP_CXX_FLAGS}>
    $<$<CXX_COMPILER_ID:Intel>:${OpenMP_CXX_FLAGS}>
  )
endif ()

#----------------------------------------------------------------------------
# Common compile options and definitions
#----------------------------------------------------------------------------
# target_compile_options(arcane_alina PUBLIC
#   # Compiler is GNU (g++):
#   $<$<CXX_COMPILER_ID:GNU>:$<BUILD_INTERFACE:-Wno-conversion;-Wno-float-conversion>>
#   # Compiler is Clang:
#   $<$<CXX_COMPILER_ID:Clang>:$<BUILD_INTERFACE:-Wno-conversion;-Wno-float-conversion>>
#   # Compiler is MSVC:
#   #$<$<CXX_COMPILER_ID:MSVC>:/bigobj>
#   #$<$<CXX_COMPILER_ID:MSVC>:/wd4715>
# )

# target_compile_definitions(amgcl INTERFACE
#   # Compiler is MSVC:
#     $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>
#     $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>
#     $<$<CXX_COMPILER_ID:MSVC>:_VARIADIC_MAX=10>
#     $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
#   )

#----------------------------------------------------------------------------
# Eigen backend
#----------------------------------------------------------------------------
find_path(EIGEN_INCLUDE Eigen/SparseCore PATH_SUFFIXES eigen3)
if (EIGEN_INCLUDE)
    add_library(eigen_target INTERFACE)
    target_include_directories(eigen_target INTERFACE ${EIGEN_INCLUDE})
    # target_compile_options(eigen_target INTERFACE
    #     $<$<CXX_COMPILER_ID:GNU>:-Wno-int-in-bool-context>
    #     $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
    #     $<$<CXX_COMPILER_ID:GNU>:-Wno-deprecated-copy>
    #     $<$<CXX_COMPILER_ID:Clang>:-Wno-deprecated-copy>
    #     $<$<CXX_COMPILER_ID:Clang>:-Wno-c++11-long-long>
    #     )
    target_compile_definitions(eigen_target INTERFACE ARCANE_ALINA_HAVE_EIGEN)
endif()

#----------------------------------------------------------------------------
# Find CUDA
#----------------------------------------------------------------------------
if (ARCCORE_HAS_CUDA)
  #find_package(CUDA)
  find_package(CUDAToolkit REQUIRED)
  if (OPENMP_FOUND)
    list(APPEND CUDA_NVCC_FLAGS -Xcompiler ${OpenMP_CXX_FLAGS})
  endif()

  #if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	#  list(APPEND CUDA_NVCC_FLAGS ${CUDA_ARCH_FLAGS} -Wno-deprecated-gpu-targets)
  #
  #  list(APPEND CUDA_NVCC_FLAGS -Xcompiler -fPIC -Xcompiler -Wno-vla)
  #endif()

  add_library(cuda_target INTERFACE)
  target_link_libraries(cuda_target INTERFACE ${CUDA_cusparse_LIBRARY})
endif()

#----------------------------------------------------------------------------
# Find MPI
#----------------------------------------------------------------------------
find_package(MPI 3.0)
if (MPI_CXX_FOUND)
  message(STATUS "FOUND MPI_CXX for Alina")
    # Need this to comply with CMP004 policy:
    string(STRIP "${MPI_CXX_LINK_FLAGS}" MPI_CXX_LINK_FLAGS)
    string(STRIP "${MPI_CXX_LIBRARIES}" MPI_CXX_LIBRARIES)

    add_library(mpi_target INTERFACE)
    target_link_libraries(mpi_target INTERFACE arcane_alina ${MPI_CXX_LIBRARIES})
    target_include_directories(mpi_target INTERFACE ${MPI_CXX_INCLUDE_PATH})
    # target_compile_options(mpi_target INTERFACE
    #     $<$<CXX_COMPILER_ID:GNU>:-Wno-long-long>
    #     $<$<CXX_COMPILER_ID:GNU>:-Wno-literal-suffix>
    #     $<$<CXX_COMPILER_ID:Clang>:-Wno-long-long>
    #     )
endif()

#----------------------------------------------------------------------------
# Find Pastix
#----------------------------------------------------------------------------
find_package(Metis  QUIET)
find_package(BLAS   QUIET)
find_library(PTHREAD_LIBRARY pthread)
find_library(BZ2_LIBRARY bz2)
find_library(HWLOC_LIBRARY hwloc)

arcane_find_package(Parmetis)
if(Parmetis_FOUND)
  message(STATUS "Found Parmetis for Alina")
  target_compile_definitions(arcane_alina PUBLIC ARCANE_ALINA_HAVE_PARMETIS ARCANE_ALINA_HAVE_PARMETIS)
  target_link_libraries(arcane_alina PUBLIC arccon::Parmetis)
endif()

option(ARCANE_ALINA_BUILD_TESTS    ON)
option(ARCANE_ALINA_BUILD_EXAMPLES ON)

option(ARCANE_ALINA_DISABLE_RARE_COMPONENTS OFF)

if(ARCANE_ALINA_DISABLE_RARE_COMPONENTS)
  add_definitions(
    -DARCANE_ALINA_RUNTIME_DISABLE_MULTICOLOR_GS
    -DARCANE_ALINA_RUNTIME_DISABLE_SPAI1
    -DARCANE_ALINA_RUNTIME_DISABLE_CHEBYSHEV
  )
endif()

#if (ARCANE_ALINA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
#endif()

add_subdirectory(examples)
add_subdirectory(tutorial)

if (ARCANE_ALINA_BUILD_EXAMPLES)

  #----------------------------------------------------------------------------
  # Build Fortran example
  #----------------------------------------------------------------------------
  option(ARCANE_ALINA_HAVE_FORTRAN "Build Fortran example" OFF)
  if (ARCANE_ALINA_HAVE_FORTRAN)
    add_subdirectory(fortran)
  endif()

  #----------------------------------------------------------------------------
  # Build Python example
  #----------------------------------------------------------------------------
  option(ARCANE_ALINA_HAVE_PYTHON "Build Python example" OFF)
  if (ARCANE_ALINA_HAVE_PYTHON)
    add_subdirectory(pyamgcl)
  endif()
endif()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
