# Sets up the build type.
function (dnaoad_setup_build_type default_build_type)

  get_property(_multiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if (_multiConfig)
    message(STATUS "Build type multi-config (build type selected at the build step)")
  else ()
    if (NOT CMAKE_BUILD_TYPE)
      message(STATUS "Build type ${default_build_type} (default single-config)")
      set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Build type" FORCE)
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING "Choose the type of build")
      set_property(
        CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo"
      )
    else ()
      message(STATUS "DNAOAD: build type: ${CMAKE_BUILD_TYPE} (manually selected single-config)")
    endif ()
  endif ()

endfunction ()

# Replaces the extension of a given file
#
# Args:
#     oldext [in]: Old extension
#     newext [in]: New extension
#     fname [in]: File name in which extension should be replaced.
#     newfname [out]: File name after extension replacement.
#
function(dnaoad_replace_extension oldext newext fname newfname)

  string(REGEX REPLACE "\\.${oldext}$" ".${newext}" _newfname ${fname})
  set(${newfname} ${_newfname} PARENT_SCOPE)

endfunction()


# Registers files for preprocessing
#
# Args:
#     preproc [in]: Preprocessor to use
#     preprocopts [in]:  Preprocessor command line arguments (but not in/out file)
#     oldext [in]: Extension of the unpreprocessed files.
#     newext [in]: Extension of the preprocessed files.
#     oldfiles [in]: List of unpreprocessed file names.
#     newfiles [out]: List of preprocessed file names.
#
function(dnaoad_preprocess preproc preprocopts oldext newext oldfiles newfiles)

  set(_newfiles)
  foreach(oldfile IN LISTS oldfiles)
    # Start with an absolulte path, so that the correct relative path is calculated thereafter
    get_filename_component(oldfile ${oldfile} ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR})
    file(RELATIVE_PATH oldfile ${CMAKE_CURRENT_SOURCE_DIR} ${oldfile})
    dnaoad_replace_extension(${oldext} ${newext} ${oldfile} newfile)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      COMMAND ${preproc} ${preprocopts} ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile} ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile})
    list(APPEND _newfiles ${CMAKE_CURRENT_BINARY_DIR}/${newfile})
  endforeach()
  set(${newfiles} ${_newfiles} PARENT_SCOPE)

endfunction()



# Build -D command line arguments for Fypp preprocessor based on current configuration
#
# Args:
#     fyppflags [inout]: Current Fypp flags on enter, with -D options extended flags on exit.
#
function (dnaoad_add_fypp_defines fyppflags)

  include(CheckFortranSourceRuns)

  set(_fyppflags "${${fyppflags}}")

  # Check for compiler int32 support
  check_fortran_source_runs(
    "program test_int32
      use iso_fortran_env, only: int32
      if (int32 < 0) then
        ! int32 is negative if not supported by compiler under F08 standard
        error stop (int32)
      else
        stop
      end if
    end program test_int32"
    HAVE_INT32_SUPPORT
  )
  if (HAVE_INT32_SUPPORT)
    message(STATUS "Fortran compiler supports iso_fortran_env int32.")
    list(APPEND _fyppflags -DINT32=1)
  else()
    list(APPEND _fyppflags -DINT32=0)
    message(WARNING "Fortran compiler does not support iso_fortran_env int32.")
  endif()

  # Check for compiler int64 support
  check_fortran_source_runs(
    "program test_int64
      use iso_fortran_env, only: int64
      if (int64 < 0) then
        ! int64 is negative if not supported by compiler under F08 standard
        error stop (int64)
      else
        stop
      end if
    end program test_int64"
    HAVE_INT64_SUPPORT
  )
  if (HAVE_INT64_SUPPORT)
    message(STATUS "Fortran compiler supports iso_fortran_env int64.")
    list(APPEND _fyppflags -DINT64=1)
  else()
    list(APPEND _fyppflags -DINT64=0)
    message(WARNING "Fortran compiler does not support iso_fortran_env int64.")
  endif()

  # Check for compiler real32 support
  check_fortran_source_runs(
    "program test_real32
      use iso_fortran_env, only: real32
      if (real32 < 0) then
        ! real32 is negative if not supported by compiler under F08 standard
        error stop (real32)
      else
        stop
      end if
    end program test_real32"
    HAVE_REAL32_SUPPORT
  )
  if (HAVE_REAL32_SUPPORT)
    message(STATUS "Fortran compiler supports iso_fortran_env real32.")
    list(APPEND _fyppflags -DREAL32=1)
  else()
    list(APPEND _fyppflags -DREAL32=0)
    message(WARNING "Fortran compiler does not support iso_fortran_env real32.")
  endif()

  # Check for compiler real64 support
  check_fortran_source_runs(
    "program test_real64
      use iso_fortran_env, only: real64
      if (real64 < 0) then
        ! real64 is negative if not supported by compiler under F08 standard
        error stop (real64)
      else
        stop
      end if
    end program test_real64"
    HAVE_REAL64_SUPPORT
  )
  if (HAVE_REAL64_SUPPORT)
    message(STATUS "Fortran compiler supports iso_fortran_env real64.")
    list(APPEND _fyppflags -DREAL64=1)
  else()
    list(APPEND _fyppflags -DREAL64=0)
    message(WARNING "Fortran compiler does not support iso_fortran_env real64.")
  endif()

  # Check for compiler real128 support
  check_fortran_source_runs(
    "program test_real128
      use iso_fortran_env, only: real128
      if (real128 < 0) then
        ! real128 is negative if not supported by compiler under F08 standard
        error stop (real128)
      else
        stop
      end if
    end program test_real128"
    HAVE_REAL128_SUPPORT
  )
  if (HAVE_REAL128_SUPPORT)
    message(STATUS "Fortran compiler supports iso_fortran_env real128.")
    list(APPEND _fyppflags -DREAL128=1)
  else()
    list(APPEND _fyppflags -DREAL128=0)
    message(WARNING "Fortran compiler does not support iso_fortran_env real128.")
  endif()

  set(${fyppflags} ${_fyppflags} PARENT_SCOPE)

endfunction()
