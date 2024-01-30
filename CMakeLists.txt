cmake_minimum_required(VERSION 3.19)

project(libomb
        LANGUAGES Fortran C)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  #set(CMAKE_INSTALL_PREFIX "/usr/local")
  #set(CMAKE_INSTALL_PREFIX "/Users/rmg3/sandbox/cmake.demos")
  set(CMAKE_INSTALL_PREFIX "/u/robert.grumbine/rgdev/mmablib")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  #set(CMAKE_INSTALL_PREFIX "/my/default" CACHE PATH "..." FORCE)
  #set(CMAKE_INSTALL_PREFIX "/Users/rmg3/sandbox/cmake.demos")
  set(CMAKE_INSTALL_PREFIX "/u/robert.grumbine/rgdev/mmablib")
endif()

#set(CMAKE_INSTALL_PREFIX "/Users/rmg3/sandbox/cmake.demos")
set(CMAKE_INSTALL_PREFIX "/u/robert.grumbine/rgdev/mmablib")

add_subdirectory( sorc )