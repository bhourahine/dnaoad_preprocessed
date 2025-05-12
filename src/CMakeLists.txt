# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

add_library(dnaoad)
add_library(DNAOAD::dnaoad ALIAS dnaoad)

set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/constants.F90
  ${curdir}/dnaoad.F90
  ${curdir}/matrixfunctions.F90
  ${curdir}/multivariablefunctions.F90
  ${curdir}/types.F90
)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp})

dnaoad_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

target_sources(
  dnaoad
  PRIVATE
  ${all-sources-f90-preproc}
)

target_include_directories(
  dnaoad PUBLIC
  $<BUILD_INTERFACE:${moduledir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${DNAOAD_INSTALL_MODULEDIR}>
)

set_target_properties(
  dnaoad
  PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    EXPORT_NAME dnaoad
    OUTPUT_NAME dnaoad
    Fortran_MODULE_DIRECTORY "${moduledir}"
)

if (DNAOAD_INSTALL)
  install(
    TARGETS dnaoad
    EXPORT DNAOADTargets
    LIBRARY
      DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      COMPONENT DNAOAD_runtime
      NAMELINK_COMPONENT DNAOAD_development
    PUBLIC_HEADER
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
      COMPONENT DNAOAD_development
  )
  install(
    DIRECTORY "${moduledir}/"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/${DNAOAD_INSTALL_MODULEDIR}"
    COMPONENT DNAOAD_development
  )
endif ()
