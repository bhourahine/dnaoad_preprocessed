# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

list(APPEND examples
  dnaoad_example1
  dnaoad_example2
  dnaoad_example3
)

foreach(example IN LISTS examples)
  add_executable(${example})
  target_sources(${example} PRIVATE ${example}.f90)
  set_target_properties(
    ${example} PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}"
  )
  target_link_libraries(${example} PRIVATE dnaoad)
endforeach()
