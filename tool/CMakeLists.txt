# include guard
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(
    FATAL_ERROR
      "CMake shouldn't be executed for 'src' subdirectory, but for the top directory of mVMC."
  )
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_library(key2lower STATIC key2lower.c)
add_executable(greenr2k greenr2k.F90)
target_link_libraries(greenr2k key2lower ${LAPACK_LIBRARIES})

install(DIRECTORY DESTINATION "$ENV{HOME}/.local/share/dvmc")
install(DIRECTORY DESTINATION "$ENV{HOME}/.local/bin")
install(TARGETS greenr2k RUNTIME DESTINATION $ENV{HOME}/.local/share/dvmc)

# Scripts
configure_file(wout2geom.sh wout2geom.sh COPYONLY)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/wout2geom.sh
  DESTINATION $ENV{HOME}/.local/share/dvmc
  PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
    WORLD_READ
    WORLD_EXECUTE)

configure_file(respack2wan90.py respack2wan90.py COPYONLY)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/respack2wan90.py
  DESTINATION $ENV{HOME}/.local/share/dvmc
  PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
    WORLD_READ
    WORLD_EXECUTE)
