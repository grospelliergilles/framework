include (srcs.cmake)

#----------------------------------------------------------------------------
# Find Boost
#----------------------------------------------------------------------------
option(Boost_USE_STATIC_LIBS "Use static versions of Boost libraries" OFF)
if (WIN32)
  set(Boost_USE_STATIC_LIBS ON)
endif ()

find_package(Boost COMPONENTS
  program_options
  serialization
  unit_test_framework
)

#----------------------------------------------------------------------------
# Builtin backend
#----------------------------------------------------------------------------
add_library(amgcl INTERFACE)
add_library(amgcl::amgcl ALIAS amgcl)
target_sources(amgcl INTERFACE ${ARCANE_ALINA_SOURCES})
  
target_compile_features(amgcl INTERFACE
  cxx_auto_type
  cxx_range_for
  cxx_rvalue_references
  cxx_right_angle_brackets
  cxx_static_assert
)

if (Boost_FOUND)
  target_include_directories(amgcl SYSTEM INTERFACE
    ${Boost_INCLUDE_DIRS}
  )
else()
  target_compile_definitions(amgcl INTERFACE AMGCL_NO_BOOST)
endif()

target_include_directories(amgcl INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

find_package(OpenMP)
if (OPENMP_FOUND)
  target_compile_options(amgcl INTERFACE ${OpenMP_CXX_FLAGS})
  target_link_libraries(amgcl INTERFACE
    $<$<CXX_COMPILER_ID:GNU>:${OpenMP_CXX_FLAGS}>
    $<$<CXX_COMPILER_ID:Clang>:${OpenMP_CXX_FLAGS}>
    $<$<CXX_COMPILER_ID:Intel>:${OpenMP_CXX_FLAGS}>
  )
endif ()

#----------------------------------------------------------------------------
# Common compile options and definitions
#----------------------------------------------------------------------------
target_compile_options(amgcl INTERFACE
  # Compiler is GNU (g++):
  $<$<CXX_COMPILER_ID:GNU>:$<BUILD_INTERFACE:-Wall;-Wextra;-Wpedantic>>
  # Compiler is Clang:
  $<$<CXX_COMPILER_ID:Clang>:$<BUILD_INTERFACE:-Wall;-Wextra;-Wpedantic;-Wno-c99-extensions>>
  # Compiler is MSVC:
  $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
  $<$<CXX_COMPILER_ID:MSVC>:/wd4715>
)

target_compile_definitions(amgcl INTERFACE
  # Compiler is MSVC:
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>
    $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>
    $<$<CXX_COMPILER_ID:MSVC>:_VARIADIC_MAX=10>
    $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
  )

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
    target_compile_definitions(eigen_target INTERFACE AMGCL_HAVE_EIGEN)
endif()

#----------------------------------------------------------------------------
# Find CUDA
#----------------------------------------------------------------------------
# if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
#     find_package(CUDA)
#     if (CUDA_FOUND)
#         set(CUDA_TARGET_ARCH "Kepler Maxwell Pascal" CACHE STRING "Target architecture(s) for CUDA")
#         cuda_select_nvcc_arch_flags(CUDA_ARCH_FLAGS ${CUDA_TARGET_ARCH})

#         if (OPENMP_FOUND)
#             list(APPEND CUDA_NVCC_FLAGS -Xcompiler ${OpenMP_CXX_FLAGS})
#         endif()

#         if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
# 	    list(APPEND CUDA_NVCC_FLAGS
#                 ${CUDA_ARCH_FLAGS} -Wno-deprecated-gpu-targets)

#             list(APPEND CUDA_NVCC_FLAGS -Xcompiler -fPIC -Xcompiler -Wno-vla)
#         endif()

#         add_library(cuda_target INTERFACE)
#         target_link_libraries(cuda_target INTERFACE ${CUDA_cusparse_LIBRARY})
#     endif()
# endif()

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
    target_link_libraries(mpi_target INTERFACE amgcl ${MPI_CXX_LIBRARIES})
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

option(AMGCL_BUILD_TESTS    ON)
option(AMGCL_BUILD_EXAMPLES ON)

option(AMGCL_DISABLE_RARE_COMPONENTS OFF)

if(AMGCL_DISABLE_RARE_COMPONENTS)
  add_definitions(
    -DAMGCL_RUNTIME_DISABLE_MULTICOLOR_GS
    -DAMGCL_RUNTIME_DISABLE_SPAI1
    -DAMGCL_RUNTIME_DISABLE_CHEBYSHEV
  )
endif()

#if (AMGCL_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
#endif()

add_subdirectory(lib)
add_subdirectory(examples)
add_subdirectory(tutorial)

if (AMGCL_BUILD_EXAMPLES)

  #----------------------------------------------------------------------------
  # Build Fortran example
  #----------------------------------------------------------------------------
  option(AMGCL_HAVE_FORTRAN "Build Fortran example" OFF)
  if (AMGCL_HAVE_FORTRAN)
    add_subdirectory(fortran)
  endif()

  #----------------------------------------------------------------------------
  # Build Python example
  #----------------------------------------------------------------------------
  option(AMGCL_HAVE_PYTHON "Build Python example" OFF)
  if (AMGCL_HAVE_PYTHON)
    add_subdirectory(pyamgcl)
  endif()
endif()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
