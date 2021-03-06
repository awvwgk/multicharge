# This file is part of multicharge.
# SPDX-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)
get_directory_property(is-subproject PARENT_DIRECTORY)

project(
  "multicharge"
  LANGUAGES "Fortran"
  VERSION "0.2.0"
  DESCRIPTION "Molecular structure store for testing"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Collect subprojects
set(lib-deps)
add_subdirectory("subprojects" EXCLUDE_FROM_ALL)

# General configuration information
add_subdirectory("config")

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

# Collect source of the project
set(srcs)
add_subdirectory("src")

# multicharge library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
  "LAPACK::LAPACK"
  "$<$<VERSION_LESS:${CMAKE_VERSION},3.20>:BLAS::BLAS>"
)
if(WITH_OpenMP)
  target_link_libraries(
    "${PROJECT_NAME}-lib"
    PUBLIC
    "OpenMP::OpenMP_Fortran"
  )
endif()
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

# Add example application
add_subdirectory("app")

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
)
# Package license files
install(
  FILES
  "LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# add the testsuite
enable_testing()
add_subdirectory("test")
