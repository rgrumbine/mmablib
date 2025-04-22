cmake_minimum_required(VERSION 3.19)

project(libomb
        LANGUAGES Fortran C)

#message("on entry HOME " HOME="${HOME}")
#message("on entry mmablib " CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  #set(CMAKE_INSTALL_PREFIX "/usr/local")
  #set(CMAKE_INSTALL_PREFIX "/Users/rmg3/sandbox/cmake.demos")
  message("CMAKE_INSTALL_PREFIX not defined")
  set(CMAKE_INSTALL_PREFIX ${HOME}/mmablib)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  #set(CMAKE_INSTALL_PREFIX "/my/default" CACHE PATH "..." FORCE)
  message("CMAKE_INSTALL_PREFIX not default")
  set(CMAKE_INSTALL_PREFIX ${HOME}/mmablib)
endif()

message("mmablib " CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
message("mmablib " CMAKE_INSTALL_LIBDIR="${CMAKE_INSTALL_LIBDIR}")

add_subdirectory( sorc )
