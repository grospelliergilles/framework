﻿# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Backend Arcane pour CUDA

set(ARCANE_SOURCES
  CudaAccelerator.cc
  CudaAccelerator.h
  )

arcane_find_package(CUDAToolkit REQUIRED)

# Créé une cible interface pour propager les options de compilation
# communes pour la compilation CUDA

add_library(arcane_cuda_compile_flags INTERFACE)

option(ARCANE_CUDA_DEVICE_DEBUG "If True, add '--device-debug' to cuda compiler flags" OFF)
set(_CUDA_DEBUG_FLAGS "-lineinfo")
if (ARCANE_CUDA_DEVICE_DEBUG)
  set(_CUDA_DEBUG_FLAGS "--device-debug")
endif()

target_compile_options(arcane_cuda_compile_flags INTERFACE
  "$<$<COMPILE_LANGUAGE:CUDA>:-ccbin>"
  "$<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CXX_COMPILER}>"
  "$<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>"
  "$<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>"
  "$<$<COMPILE_LANGUAGE:CUDA>:${_CUDA_DEBUG_FLAGS}>"
  "$<$<COMPILE_LANGUAGE:CUDA>:-g>"
  "$<$<COMPILE_LANGUAGE:CUDA>:-Werror>"
  "$<$<COMPILE_LANGUAGE:CUDA>:cross-execution-space-call>"
)

install(TARGETS arcane_cuda_compile_flags EXPORT ArcaneTargets)

arcane_add_library(arcane_accelerator_cuda
  INPUT_PATH ${Arcane_SOURCE_DIR}/src
  RELATIVE_PATH arcane/accelerator/cuda
  FILES ${ARCANE_SOURCES}
  )

target_link_libraries(arcane_accelerator_cuda PUBLIC
  arcane_core arcane_cuda_compile_flags CUDA::cudart
  )

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Implémentation des routines spécifiques CUDA.
# Pour que cela fonctionne, cette bibliothèque doit être statique
# (sinon, il semble que CMake ne l'ajoute pas à l'édition de lien des cibles
# qui en dépende mais peut-être s'agit-il d'une limitation de 'nvcc')
# Elle permet aussi de propager les options de compilation aux utilisateurs
# de cette bibliothèque
add_library(arcane_accelerator_cuda_impl STATIC
  Reduce.cu
#  CudaReduceImpl.h
  )
set_target_properties(arcane_accelerator_cuda_impl PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(arcane_accelerator_cuda_impl PROPERTIES LINKER_LANGUAGE CUDA)
target_link_libraries(arcane_accelerator_cuda_impl PUBLIC arcane_core arcane_cuda_compile_flags)
target_link_libraries(arcane_accelerator_cuda_impl PUBLIC CUDA::cudart)

# ----------------------------------------------------------------------------

arcane_register_library(arcane_accelerator_cuda  OPTIONAL)

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
