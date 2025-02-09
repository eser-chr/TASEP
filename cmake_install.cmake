# Install script for directory: /home/chris/Desktop/TASEP

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/chris/anaconda3/lib/python3.12/site-packages" TYPE MODULE FILES "/home/chris/Desktop/TASEP/tasep.so")
  if(EXISTS "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/chris/anaconda3/lib/python3.12/site-packages" TYPE MODULE FILES "/home/chris/Desktop/TASEP/neighbors.so")
  if(EXISTS "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/neighbors.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/chris/anaconda3/lib/python3.12/site-packages" TYPE MODULE FILES "/home/chris/Desktop/TASEP/tasep_profile.so")
  if(EXISTS "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/chris/anaconda3/lib/python3.12/site-packages/tasep_profile.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/chris/Desktop/TASEP/external/bucket/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
  file(WRITE "/home/chris/Desktop/TASEP/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
