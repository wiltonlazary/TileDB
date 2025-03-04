#
# test/regression/CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2022 TileDB, Inc.
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

set(TILEDB_CORE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
list(APPEND TILEDB_CORE_INCLUDE_DIR "${CMAKE_BINARY_DIR}/tiledb/")

find_package(Catch_EP REQUIRED)

set(SOURCES
  targets/sc-12024.cc
  targets/sc-15387.cc
  targets/sc-18836.cc
)

if (TILEDB_SERIALIZATION)
  list(APPEND SOURCES targets/sc-18250.cc)
endif()

if (TILEDB_CPP_API)
  list(APPEND SOURCES targets/sc-19240_cppapi-vfs-exception.cc)
endif()

add_executable(tiledb_regression
  EXCLUDE_FROM_ALL
  regression.cc
  ${SOURCES}
)

if (NOT MSVC)
  target_compile_options(tiledb_regression PRIVATE -Wno-deprecated-declarations)
endif()

target_link_libraries(tiledb_regression
  PUBLIC
    Catch2::Catch2WithMain
    local_install
    tiledb_shared
)

target_include_directories(tiledb_regression
  PRIVATE
    ${CMAKE_INSTALL_PREFIX}/include
    "${TILEDB_CORE_INCLUDE_DIR}"
    ${TILEDB_EP_INSTALL_PREFIX}/include
)

add_test(
  NAME "tiledb_regression"
  COMMAND $<TARGET_FILE:tiledb_regression> --durations=yes
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
