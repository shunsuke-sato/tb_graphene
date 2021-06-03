set(EXECUTABLE tb_gra)

add_executable(${EXECUTABLE}
  main.f90
  parallel.f90
  communication.f90
  inputoutput.f90
  )



set_target_properties(${EXECUTABLE}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )

