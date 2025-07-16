# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

list(APPEND tests
  testarithmetic
  testexponentiation
  testjacobian
  testmatrix
  testtrig
  testtype
  testvector
)

foreach(test IN LISTS tests)
  add_executable(${test})
  target_sources(${test} PRIVATE ${test}.f90)
  set_target_properties(
    ${test} PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}"
  )
  target_link_libraries(${test} PRIVATE dnaoad)
  target_link_libraries(${test} PRIVATE Fortuno::fortuno_serial)
  
  add_test(
    NAME ${test}
    COMMAND ${test}
  )
endforeach()
