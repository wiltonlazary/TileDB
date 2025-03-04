#
# CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2017-2021 TileDB, Inc.
# Copyright (c) 2016 MIT and Intel Corporation
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

############################################################
# CMake setup
############################################################

cmake_minimum_required(VERSION 3.21)

############################################################
# Parse version file
# credit: https://stackoverflow.com/a/47084079

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/tiledb/sm/c_api/tiledb_version.h" VERFILE)
if (NOT VERFILE)
  message(FATAL_ERROR "Failed to parse tiledb_version.h!")
endif()

string(REGEX MATCH "TILEDB_VERSION_MAJOR ([0-9])*" _ ${VERFILE})
set(TILEDB_VERSION_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "TILEDB_VERSION_MINOR ([0-9]+)*" _ ${VERFILE})
set(TILEDB_VERSION_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "TILEDB_VERSION_PATCH ([0-9]+)*" _ ${VERFILE})
set(TILEDB_VERSION_PATCH ${CMAKE_MATCH_1})

set(TILEDB_VERSION "${TILEDB_VERSION_MAJOR}.${TILEDB_VERSION_MINOR}.${TILEDB_VERSION_PATCH}")
############################################################
# Check for regex characters in the most important paths
# fixes https://github.com/TileDB-Inc/TileDB/issues/1799
option(TILEDB_ALLOW_REGEX_CHAR_PATH "If true, allow regex characters in source, build, or install path." FALSE)
mark_as_advanced(TILEDB_ALLOW_REGEX_CHAR_PATH)
set(REGEX_CHARS "[\\^\\$\\+\\*\\?\\|\\(\\)]") # note: must be escaped, and regex doesn't work with \[\] entries
set(REGEX_CHAR_PATH_MSG " contains a REGEX character and may break CMakeList processing. Please use"
                        " a different path, or set TILEDB_ALLOW_REGEX_CHAR_PATH to override.")
if (NOT TILEDB_ALLOW_REGEX_CHAR_PATH)
  if (CMAKE_CURRENT_SOURCE_DIR MATCHES ${REGEX_CHARS})
    message(FATAL_ERROR "CMAKE_CURRENT_SOURCE_DIR ${REGEX_CHAR_PATH_MSG}:\n  '${CMAKE_CURRENT_SOURCE_DIR}'")
  elseif (CMAKE_CURRENT_SOURCE_DIR MATCHES ${REGEX_CHARS})
    message(FATAL_ERROR "CMAKE_CURRENT_BINARY_DIR ${REGEX_CHAR_PATH_MSG}:\n  '${CMAKE_CURRENT_BINARY_DIR}'")
  elseif (CMAKE_CURRENT_SOURCE_DIR MATCHES ${REGEX_CHARS})
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX ${REGEX_CHAR_PATH_MSG}:\n  '${CMAKE_INSTALL_PREFIX}'")
  endif()
endif()

############################################################
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Options")

set(TILEDB_CMAKE_INPUTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/inputs")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(APPLE)
  # Use @rpath on macOS for building shared libraries.
  set(CMAKE_MACOSX_RPATH ON)
  # Don't allow macOS .frameworks to be used for dependencies.
  set(CMAKE_FIND_FRAMEWORK NEVER)
endif()

# Set C++17 as required standard for all C++ targets.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
  # Use GNU extensions under Cygwin
  set(CMAKE_CXX_EXTENSIONS ON)
else()
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# Set -fvisibility=hidden (or equivalent) flags by default.
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

############################################################
# Build options
############################################################

# Note: when adding options, make sure to forward them via INHERITED_CMAKE_ARGS
# in TileDB-Superbuild.cmake.

option(TILEDB_SUPERBUILD "If true, perform a superbuild (builds all missing dependencies)." ON)
option(TILEDB_FORCE_ALL_DEPS "If true, force superbuild to download and build all dependencies, even those installed on the system." OFF)
option(TILEDB_VERBOSE "Prints TileDB errors with verbosity" OFF)
option(TILEDB_S3 "Enables S3/minio support using aws-cpp-sdk" OFF)
option(TILEDB_AZURE "Enables Azure Storage support using azure-storage-cpp" OFF)
option(TILEDB_GCS "Enables GCS Storage support using google-cloud-cpp" OFF)
option(TILEDB_HDFS "Enables HDFS support using the official Hadoop JNI bindings" OFF)
option(TILEDB_WERROR "Enables the -Werror flag during compilation." ON)
option(TILEDB_ASSERTIONS "Build with assertions enabled (default off for release, on for debug build)." OFF)
option(TILEDB_CPP_API "Enables building of the TileDB C++ API." ON)
option(TILEDB_CMAKE_IDE "(Used for CLion builds). Disables superbuild and sets the EP install dir." OFF)
option(TILEDB_STATS "Enables internal TileDB statistics gathering." ON)
option(TILEDB_STATIC "Enables building TileDB as a static library." OFF)
option(TILEDB_TESTS "If true, enables building the TileDB unit test suite" ON)
option(TILEDB_TOOLS "If true, enables building the TileDB tools" OFF)
option(TILEDB_SERIALIZATION "If true, enables building with support for query serialization" OFF)
option(TILEDB_CCACHE "If true, enables use of 'ccache' (if present)" OFF)
option(TILEDB_ARROW_TESTS "If true, enables building the arrow adapter unit tests" OFF)
option(TILEDB_CRC32 "If true, enables building crc32 and a simple linkage test" OFF)
option(TILEDB_WEBP "If true, enables building webp and a simple linkage test" OFF)
option(TILEDB_ABSEIL "If true, enables building abseil and a simple linkage test" OFF)
option(TILEDB_LOG_OUTPUT_ON_FAILURE "If true, print error logs if dependency sub-project build fails" ON)
option(TILEDB_SKIP_S3AWSSDK_DIR_LENGTH_CHECK "If true, skip check needed path length for awssdk (TILEDB_S3) dependent builds" OFF)
option(TILEDB_EXPERIMENTAL_FEATURES "If true, build and include experimental features" OFF)

set(TILEDB_INSTALL_LIBDIR "" CACHE STRING "If non-empty, install TileDB library to this directory instead of CMAKE_INSTALL_LIBDIR.")

# early WIN32 audit of path length for aws sdk build where
# insufficient available path length causes sdk build failure.
if (WIN32 AND NOT TILEDB_SKIP_S3AWSSDK_DIR_LENGTH_CHECK)
  if (TILEDB_SUPERBUILD)
    if (TILEDB_S3)
      string(LENGTH ${CMAKE_CURRENT_BINARY_DIR} LENGTH_CMAKE_CURRENT_BINARY_DIR)
      if ( NOT (LENGTH_CMAKE_CURRENT_BINARY_DIR LESS 61))
        message(FATAL_ERROR " build directory path likely too long for building awssdk/dependencies!")
      return()
      endif()
    endif()
  endif()
endif()

# enable assertions by default for debug builds
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(TILEDB_ASSERTIONS TRUE)
endif()
include(TileDBAssertions)

############################################################
# Superbuild setup
############################################################

# Set the variable used when calling find_package(), find_file() etc.
# to determine if NO_DEFAULT_PATH should be passed.
if (TILEDB_FORCE_ALL_DEPS)
  set(TILEDB_DEPS_NO_DEFAULT_PATH NO_DEFAULT_PATH)
else()
  set(TILEDB_DEPS_NO_DEFAULT_PATH)
endif()

# If this is an in-IDE build, we need to disable the superbuild and explicitly
# set the EP base dir. The normal 'cmake && make' process won't need this step,
# it is for better CLion support of this superbuild architecture.
if (TILEDB_CMAKE_IDE)
  set(TILEDB_SUPERBUILD OFF)
  set(TILEDB_EP_BASE "${CMAKE_CURRENT_BINARY_DIR}/externals")
endif()

# When building static TileDB, we also need to install any static dependencies
# built as external projects.
set(TILEDB_INSTALL_STATIC_DEPS ${TILEDB_STATIC})

# Set main TileDB Source
set(TILEDB_BASE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/tiledb")
# Set experimental TileDB Source
set(TILEDB_EXPERIMENTAL_BASE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/experimental")

# Perform superbuild config and exit.
if (TILEDB_SUPERBUILD)
  project(TileDB-Superbuild)
  message(STATUS "Starting TileDB superbuild.")
  include("cmake/TileDB-Superbuild.cmake")
  message(STATUS "Install prefix is ${CMAKE_INSTALL_PREFIX}")
  # The superbuild file incorporates this file as an external project with the
  # superbuild off. Thus we stop processing this file here, knowing that we'll
  # come back later.
  return()
endif()

project(TileDB)
message(STATUS "Starting TileDB regular build.")
# Paths to locate the installed external projects.
set(TILEDB_EP_SOURCE_DIR "${TILEDB_EP_BASE}/src")
set(TILEDB_EP_INSTALL_PREFIX "${TILEDB_EP_BASE}/install")

############################################################
# Compile options/definitions for all targets
############################################################

# Set compiler flags
if (MSVC)
  # Turn on standards-conformance mode
  add_compile_options("/permissive-")
  # We disable some warnings that are not present in gcc/clang -Wall:
  #   C4101: unreferenced local variable
  #   C4146: unary minus operator applied to unsigned type
  #   C4244: conversion warning of floating point to integer type.
  #   C4251: C++ export warning
  #   C4456: local variable hiding previous local variable
  #   C4457: local variable hiding function parameter
  #   C4702: unreachable code
  #   C4800: warning implicit cast int to bool
  add_compile_options(/W4 /wd4101 /wd4146 /wd4244 /wd4251 /wd4456 /wd4457 /wd4702 /wd4800)
  # Warnings as errors:
  if (TILEDB_WERROR)
    add_compile_options(/WX)
  endif()
  # Turn off MSVC deprecation of certain standard library functions. This allows
  # other deprecations to remain visible.
  add_compile_definitions("_CRT_SECURE_NO_WARNINGS")
  # We currently need to suppress warnings about deprecation (C4996) for two cases:
  #   1. C++ API functions that call deprecated C API functions
  #   2. two warnings in `test/src/helpers.cc` that call deprecated C API functions
  add_compile_options(/wd4996)
  # Disable GDI (which we don't need, and causes some macro
  # re-definition issues if wingdi.h is included)
  add_compile_options(/DNOGDI)
  # Add /MPn flag from CMake invocation (if defined).
  add_compile_options(${MSVC_MP_FLAG})
  # Build-specific flags
  add_compile_options(
                      "$<$<CONFIG:Debug>:/DDEBUG /Od /Zi /bigobj>"
                      "$<$<CONFIG:Release>:/DNDEBUG /Ox>"
                      "$<$<CONFIG:RelWithDebInfo>:/DNDEBUG /Ox /Zi>"
                      )
else()
  add_compile_options(-Wall -Wextra)
  if (TILEDB_WERROR)
    add_compile_options(-Werror)
  endif()
  # Build-specific flags
  if (CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-DDEBUG -O0 -g3 -ggdb3 -gdwarf-3)
  elseif (CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(-O3)
  elseif (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    add_compile_options(-DNDEBUG -O3 -g3 -ggdb3 -gdwarf-3)
  elseif (CMAKE_BUILD_TYPE MATCHES "Coverage")
    add_compile_options(-DDEBUG -g3 -gdwarf-3 --coverage)
  endif()

  # Use -Wno-literal-suffix on Linux with C++ sources.
  if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-literal-suffix>)
  endif()
endif()

# Definitions for all targets
add_definitions(-D_FILE_OFFSET_BITS=64)

# AVX2 flag
include(CheckAVX2Support)
CheckAVX2Support()
if (COMPILER_SUPPORTS_AVX2)
  add_compile_options(${COMPILER_AVX2_FLAG})
endif()

#######################################################
# Header Files
#######################################################
#
# Legacy headers: All headers for both C and C++ API in a single directory.
#
# Public header structure. Users of the library see this view. External source
# code uses these names; they should be considered unchangeable.
# - `<root>/`
#     `tiledb/`
#       - C API headers
#         `tiledb.h`
#         `tiledb_experimental.h`
#          other `tiledb_*.h`
#       - C++ API headers
#         `tiledb`
#         `tiledb_experimental`
#          other `tiledb_*`
#
# Headers assume `-I<root>`. To wit, for the C API use this:
# ```C
# #include <tiledb/tiledb.h>
# ```
# For the C++ API use this:
# ```C++
# #include <tiledb/tiledb>
# ```

# Private header structure. The compiler of a user program sees this view. Only
# TileDB source uses these names. They can be changed at will, but installation
# will also have to changed, as may other aspects of the build.
# - `<root>/`
#     `tiledb/`
#         `api/`
#             `c_api/`
#                 `api_external_common.h`
#                 `tiledb_export.h` (auto-generated)
#                 `<section>/` (multiple)
#                     `<section>_api_external.h`

##################################
# C API include files
##################################
#
# FILENAME headers are copied only as their file name to ${INCLUDE_BASE}/tiledb
# RELATIVE headers are copied as their path from RELATIVE_HEADER_BASE to ${INCLUDE_BASE}
# In all cases the included path begins with `tiledb/`
# The export header is handled separatedly below.

list(APPEND TILEDB_C_API_FILENAME_HEADERS
    "${CMAKE_SOURCE_DIR}/tiledb/sm/c_api/tiledb.h"
    "${CMAKE_SOURCE_DIR}/tiledb/sm/c_api/tiledb_enum.h"
    "${CMAKE_SOURCE_DIR}/tiledb/sm/c_api/tiledb_version.h"
    "${CMAKE_SOURCE_DIR}/tiledb/sm/c_api/tiledb_experimental.h"
)
list(APPEND TILEDB_C_API_RELATIVE_HEADERS
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/api_external_common.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/config/config_api_external.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/context/context_api_external.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/error/error_api_external.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/filesystem/filesystem_api_external.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/filesystem/filesystem_api_enum.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/filter/filter_api_enum.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/filter/filter_api_external.h"
    "${CMAKE_SOURCE_DIR}/tiledb/api/c_api/filter_list/filter_list_api_external.h"
)
set(TILEDB_C_API_RELATIVE_HEADER_BASE "${CMAKE_CURRENT_SOURCE_DIR}")

##################################
# Export Header
##################################
# The export header is automatically generated later in the build (using the
# main library target). When installed, it's in the same directory as the
# compiler-visible, interface-hidden header `api_external_common.h`.

set(TILEDB_EXPORT_HEADER_DIR "${CMAKE_BINARY_DIR}/tiledb")
set(TILEDB_EXPORT_HEADER_NAME "tiledb_export.h")
set(TILEDB_EXPORT_HEADER "${TILEDB_EXPORT_HEADER_DIR}/${TILEDB_EXPORT_HEADER_NAME}")
set(TILEDB_EXPORT_HEADER_LOCALINSTALL_PATH "tiledb/api/c_api/${TILEDB_EXPORT_HEADER_NAME}")

##################################
# Local install
##################################
#
# User code within this project requires manual installation that mimics a
# an installation from a distribution artifact. It's required by examples,
# integration tests, and regression tests.
#
# Caveat: The "installation" is done at _configuration_ time, not at _build_ time.
# This is a side-effect of using `configure_file`. It's not necessary and it may
# be desirable to change it to operate at build time in the future.

set(TILEDB_LOCALINSTALL_DIR "${CMAKE_BINARY_DIR}/dist-in-build")
set(TILEDB_LOCALINSTALL_INCLUDE "${TILEDB_LOCALINSTALL_DIR}/include")

#####################
# C API
#####################
foreach(HEADER ${TILEDB_C_API_FILENAME_HEADERS})
  cmake_path(GET HEADER FILENAME HEADER_STRIPPED)
  configure_file(${HEADER} ${TILEDB_LOCALINSTALL_INCLUDE}/tiledb/${HEADER_STRIPPED} COPYONLY)
endforeach()
foreach(HEADER ${TILEDB_C_API_RELATIVE_HEADERS})
  cmake_path(RELATIVE_PATH HEADER
      BASE_DIRECTORY ${TILEDB_C_API_RELATIVE_HEADER_BASE}
      OUTPUT_VARIABLE HEADER_STRIPPED
      )
  configure_file(${HEADER} ${TILEDB_LOCALINSTALL_INCLUDE}/${HEADER_STRIPPED} COPYONLY)
endforeach()
# `configure_file` for the export header happens after it's been generated


#####################
# C++ API
#####################
file(GLOB TILEDB_CPP_HEADERS
    "${CMAKE_SOURCE_DIR}/tiledb/sm/cpp_api/*.h"
    "${CMAKE_SOURCE_DIR}/tiledb/sm/cpp_api/tiledb"
    "${CMAKE_SOURCE_DIR}/tiledb/sm/cpp_api/tiledb_experimental"
    )
foreach(HEADER ${TILEDB_CPP_HEADERS})
  string(REGEX
      REPLACE "^${CMAKE_SOURCE_DIR}/tiledb/sm/cpp_api/" ""
      HEADER_STRIPPED ${HEADER}
      )
  configure_file(${HEADER} ${TILEDB_LOCALINSTALL_INCLUDE}/tiledb/${HEADER_STRIPPED} COPYONLY)
endforeach()

#####################
# Interface library
#####################
add_library(local_install INTERFACE)
target_include_directories(local_install INTERFACE ${TILEDB_LOCALINSTALL_INCLUDE})

############################################################
# Enable testing and add subdirectories
############################################################

# Enable testing
enable_testing()

# Aggregate all link-completeness targets
add_custom_target(all_link_complete)

# Build the TileDB library experimental features
add_subdirectory(experimental)

# Build the TileDB library
add_subdirectory(tiledb)

# Build examples
add_subdirectory(examples)

# Build unit tests
if (TILEDB_TESTS)
  add_subdirectory(test/support)
  add_subdirectory(test)

  # Add cmake target for "tests" to build all unit tests executables
  add_custom_target(tests)
  add_dependencies(tests tiledb_unit)
  add_dependencies(tests test_assert)

  add_dependencies(tests unit_buffer unit_datum unit_dynamic_memory)
  add_dependencies(tests unit_exception unit_interval unit_thread_pool unit_experimental)
  add_dependencies(tests unit_array_schema unit_filter_create unit_filter_pipeline unit_metadata)
  add_dependencies(tests unit_compressors unit_query unit_misc unit_vfs unit_array)
  add_dependencies(tests unit_config)
  add_dependencies(tests unit_range_subset)
  add_dependencies(tests unit_add_ranges_list)
  add_dependencies(tests unit_consistency)
  add_dependencies(tests unit_range)
  add_dependencies(tests unit_rtree)

  # C API support
  add_dependencies(tests unit_capi_handle unit_capi_exception_wrapper)
  # C API basics
  add_dependencies(tests unit_capi_error unit_capi_config unit_capi_context)
  # C API array schema
  add_dependencies(tests unit_capi_filter unit_capi_filter_list)

  add_dependencies(tests unit_delete_update_condition)
  add_dependencies(tests unit_ast)
  add_dependencies(tests unit_serializers)

  if (TILEDB_CRC32)
    add_dependencies(tests unit_link_crc32)
  endif()

  if (TILEDB_WEBP)
    add_dependencies(tests unit_link_webp)
  endif()

  if (TILEDB_ABSEIL)
    add_dependencies(tests unit_link_absl)
  endif()

  if (ENABLE_MAGIC_TEST)
    add_dependencies(tests unit_mgc_dict)
  endif()

  add_subdirectory(test/regression)
endif()

# Build tools
if (TILEDB_TOOLS)
  add_subdirectory(tools)
endif()

###########################################################
# Uninstall
###########################################################

set(CMD "xargs printf -- '-- Uninstalling: %s\\\\n' <install_manifest.txt")
add_custom_target(
   uninstall
   COMMAND echo "Uninstalling TileDB from ${CMAKE_INSTALL_PREFIX}..."
   COMMAND eval "${CMD}"
   COMMAND xargs rm -f < install_manifest.txt
   COMMAND rmdir "${CMAKE_INSTALL_PREFIX}/include/tiledb"
   COMMAND echo "TileDB uninstalled"
)
