# Include guard
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(
    FATAL_ERROR
      "CMake shouldn't be executed for 'src' subdirectory, but for the top directory of mVMC."
  )
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

install(
  CODE "execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/dvmc $ENV{HOME}/.local/bin/dvmc
)")

message("-- dvmc symbolic link installed in: $HOME/.local/bin/dvmc")
