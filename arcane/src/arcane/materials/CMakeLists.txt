message(STATUS "Parsing ${CMAKE_CURRENT_SOURCE_DIR}")

include(srcs.cmake)

arcane_add_library(arcane_materials
  INPUT_PATH ${Arcane_SOURCE_DIR}/src
  RELATIVE_PATH arcane/materials
  FILES ${ARCANE_SOURCES}
  )

target_compile_definitions(arcane_materials PRIVATE ARCANE_COMPONENT_arcane_materials)
target_include_directories(arcane_materials PUBLIC $<BUILD_INTERFACE:${Arcane_SOURCE_DIR}/src> $<INSTALL_INTERFACE:include>)

target_link_libraries(arcane_materials PUBLIC arcane_core)
arcane_register_library(arcane_materials)
