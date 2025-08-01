# Variables influencing how subproject is obtained
set(CMAKE_REQUIRE_FIND_PACKAGE_Fypp ${DNAOAD_SUBPROJECT_REQUIRE_FIND})
set(CMAKE_DISABLE_FIND_PACKAGE_Fypp ${DNAOAD_SUBPROJECT_DISABLE_FIND})
# set FETCHCONTENT_SOURCE_DIR_FYPP to use a local copy of the subproject

# Make subproject available
FetchContent_Declare(
  Fypp
  GIT_REPOSITORY "https://github.com/aradi/fypp.git"
  GIT_TAG "main"
)
FetchContent_MakeAvailable(Fypp)

if (Fypp_FOUND)
  message(STATUS "Subproject Fypp: using installed version")
  set(FYPP "fypp")
else ()
  message(STATUS "Subproject Fypp: using local copy in ${fypp_BINARY_DIR}")
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  if(WIN32)
    if(EXISTS "${fypp_SOURCE_DIR}/fypp")
      set(FYPP "${Python3_EXECUTABLE}" "${fypp_SOURCE_DIR}/fypp")
    else ()
      set(FYPP "${Python3_EXECUTABLE}" "${fypp_SOURCE_DIR}/bin/fypp")
    endif()
  else ()
    if(EXISTS "${fypp_SOURCE_DIR}/fypp")
      set(FYPP "${fypp_SOURCE_DIR}/fypp")
    else ()
      set(FYPP "${fypp_SOURCE_DIR}/bin/fypp")
    endif()
  endif()
endif ()
