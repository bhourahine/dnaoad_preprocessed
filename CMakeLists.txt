#[=============================================================================[
#                       Basic project definition                              #
]=============================================================================]

cmake_minimum_required(VERSION 3.24...3.29)

list(APPEND CMAKE_MESSAGE_CONTEXT DNAOAD)

project(
  DNAOAD
  VERSION 0.1.0
  DESCRIPTION "Dual Numbers for Arbitrary Order Automatic Differentiation"
  LANGUAGES Fortran
)

#[=============================================================================[
#                       Options                                               #
]=============================================================================]

option(
  DNAOAD_BUILD_SHARED_LIBS
  "DNAOAD: Build as shared library"
  OFF
)

option(
  DNAOAD_BUILD_EXAMPLES
  "DNAOAD: Whether to build the examples"
  ${PROJECT_IS_TOP_LEVEL}
)
option(
  DNAOAD_BUILD_TESTS
  "DNAOAD: Whether to build the tests"
  ${PROJECT_IS_TOP_LEVEL}
)
option(
  DNAOAD_INSTALL
  "DNAOAD: Install project"
  ${PROJECT_IS_TOP_LEVEL}
)
set(
  DNAOAD_INSTALL_MODULEDIR
  "modules"
  CACHE STRING
  "DNAOAD: Sub-directory for installed Fortran module files (relative to CMAKE_INSTALL_LIBDIR)"
)

option(
  DNAOAD_SUBPROJECT_REQUIRE_FIND
  "DNAOAD: Require find_package for all subprojects"
  OFF
)
option(
  DNAOAD_SUBPROJECT_DISABLE_FIND
  "DNAOAD: Disable find_package for all subprojects"
  OFF
)

#[=============================================================================[
#                  Project configuration                                      #
]=============================================================================]

include(FetchContent)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(DNAOADHelpers)

set(BUILD_SHARED_LIBS ${DNAOAD_BUILD_SHARED_LIBS})
dnaoad_setup_build_type("RelWithDebInfo")

include(subprojects/fypp.cmake)

# Extra flags from cmake options and compiler capabilities
dnaoad_add_fypp_defines(FYPP_FLAGS)

set(FYPP_FLAGS "${FYPP_FLAGS}")

set(FYPP_CONFIG_FLAGS "${FYPP_FLAGS}")
# Make sure, the line-marker option is not set
list(REMOVE_ITEM FYPP_CONFIG_FLAGS "-n")
set(FYPP_BUILD_FLAGS "${FYPP_FLAGS}" "--file-var-root=${CMAKE_SOURCE_DIR}"
  "$<IF:$<CONFIG:Debug>,-DDEBUG=1,-DDEBUG=0>")

set(PYTHON_INTERPRETER "python3" CACHE STRING
  "Python interpreter to use for preprocessor")

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${VERSION}'")

#[=============================================================================[
#                       Main definition                                       #
]=============================================================================]

add_subdirectory(src)

if (DNAOAD_BUILD_EXAMPLES)
  add_subdirectory(example)
endif ()
if (DNAOAD_BUILD_TESTS)
  include(subprojects/Fortuno.cmake)
  enable_testing()
  add_subdirectory(test)
endif ()


#[=============================================================================[
#                       Install or Export                                     #
]=============================================================================]

if (DNAOAD_INSTALL)

  # pkg-config files
  configure_file(cmake/dnaoad.pc.in dnaoad.pc @ONLY)
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/dnaoad.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    COMPONENT DNAOAD_development
  )

  # cmake export files
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DNAOADConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )
  configure_package_config_file(
      cmake/DNAOADConfig.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/DNAOADConfig.cmake
      INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DNAOAD
  )
  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/DNAOADConfigVersion.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/DNAOADConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DNAOAD
    COMPONENT DNAOAD_development
  )

  export(
    EXPORT DNAOADTargets
    FILE DNAOADTargets.cmake
    NAMESPACE DNAOAD::
  )
  install(
    EXPORT DNAOADTargets
    FILE DNAOADTargets.cmake
    NAMESPACE DNAOAD::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DNAOAD
    COMPONENT DNAOAD_development
  )
endif ()

# Make project available for FetchContent
if (NOT PROJECT_IS_TOP_LEVEL)
  # Propagate variables
  set(DNAOAD_VERSION ${DNAOAD_VERSION} PARENT_SCOPE)
  set(DNAOAD_VERSION_MAJOR ${DNAOAD_VERSION_MAJOR} PARENT_SCOPE)
  set(DNAOAD_VERSION_MINOR ${DNAOAD_VERSION_MINOR} PARENT_SCOPE)
  set(DNAOAD_VERSION_PATCH ${DNAOAD_VERSION_PATCH} PARENT_SCOPE)
  set(DNAOAD_VERSION_TWEAK ${DNAOAD_VERSION_TWEAK} PARENT_SCOPE)
endif ()
