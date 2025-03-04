#
# experimental/tiledb/common/ports/CMakeLists.txt
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

include(common NO_POLICY_SCOPE)

option(CMAKE_EXPORT_COMPILE_COMMANDS "cmake compile commands" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CMAKE_EXPORT_COMPILE_COMMANDS "cmake compile commands" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND PORT_SOURCES
    ports.cc
)
gather_sources(${PORT_SOURCES})

#
# Object library for other units to depend upon
#
add_library(ports OBJECT ${PORT_SOURCES})
target_link_libraries(ports PUBLIC baseline $<TARGET_OBJECTS:baseline>)

#
# Test-compile of object library ensures link-completeness
#
add_executable(compile_ports EXCLUDE_FROM_ALL)
add_dependencies(all_link_complete compile_ports)
target_link_libraries(compile_ports PRIVATE ports)
target_sources(compile_ports PRIVATE test/compile_ports_main.cc)

if (TILEDB_TESTS)
    add_executable(unit_ports EXCLUDE_FROM_ALL)
    target_link_libraries(unit_ports PUBLIC ports)
    find_package(Catch_EP REQUIRED)
    target_link_libraries(unit_ports PUBLIC Catch2::Catch2)
    target_link_libraries(unit_ports PUBLIC $<TARGET_OBJECTS:thread_pool>)

    # Sources for code elsewhere required for tests
    target_sources(unit_ports PUBLIC ${DEPENDENT_SOURCES})

    # Sources for tests
    target_sources(unit_ports PUBLIC
            test/main.cc
            test/unit_ports.cc
            )

    add_test(
            NAME "unit_ports"
            COMMAND $<TARGET_FILE:unit_ports> --durations=yes
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

endif()



list(APPEND PSEUDO_NODES_SOURCES
    test/pseudo_nodes.cc
)
gather_sources(${PSEUDO_NODES_SOURCES})

#
# Object library for other units to depend upon
#
#add_library(pseudo_nodes OBJECT ${PSEUDO_NODES_SOURCES})
#target_link_libraries(pseudo_nodes PUBLIC baseline $<TARGET_OBJECTS:baseline>)

#
# Test-compile of object library ensures link-completeness
#
add_executable(compile_pseudo_nodes EXCLUDE_FROM_ALL)
target_link_libraries(compile_pseudo_nodes PRIVATE pseudo_nodes)
target_sources(compile_pseudo_nodes PRIVATE test/compile_pseudo_nodes_main.cc)

if (TILEDB_TESTS)
    add_executable(unit_pseudo_nodes EXCLUDE_FROM_ALL)
#    target_link_libraries(unit_pseudo_nodes PUBLIC pseudo_nodes)
    find_package(Catch_EP REQUIRED)
    target_link_libraries(unit_pseudo_nodes PUBLIC Catch2::Catch2)
#    target_link_libraries(unit_pseudo_nodes PUBLIC $<TARGET_OBJECTS:thread_pool>)

    # Sources for code elsewhere required for tests
    target_sources(unit_pseudo_nodes PUBLIC ${DEPENDENT_SOURCES})

    # Sources for tests
    target_sources(unit_pseudo_nodes PUBLIC
            test/main.cc
            test/unit_pseudo_nodes.cc
            )

    add_test(
            NAME "unit_pseudo_nodes"
            COMMAND $<TARGET_FILE:unit_pseudo_nodes> --durations=yes
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()


# list(APPEND SIEVE_SOURCES
#    test/sieve.cc
# )
# gather_sources(${SIEVE_SOURCES})

#
# Object library for other units to depend upon
#
# add_library(sieve OBJECT ${SIEVE_SOURCES})
# target_link_libraries(sieve PUBLIC baseline $<TARGET_OBJECTS:baseline>)

#
# Since there are no components, there is no object library to test compile
#
# add_executable(compile_sieve EXCLUDE_FROM_ALL)
# target_link_libraries(compile_sieve PRIVATE sieve)
# target_sources(compile_sieve PRIVATE test/compile_sieve_main.cc)

if (TILEDB_TESTS)
    add_executable(sieve EXCLUDE_FROM_ALL)

#    target_link_libraries(sieve PUBLIC sieve)
#    find_package(Catch_EP REQUIRED)
#    target_link_libraries(sieve PUBLIC Catch2::Catch2)
#    target_link_libraries(sieve PUBLIC $<TARGET_OBJECTS:thread_pool>)

    # Sources for code elsewhere required for tests
    target_sources(sieve PUBLIC ${DEPENDENT_SOURCES})

    # Sources for tests
    target_sources(sieve PUBLIC
            test/sieve.cc
            )

# Not a test (yet) per se
#    add_test(
#            NAME "sieve"
#            COMMAND $<TARGET_FILE:sieve> --durations=yes
#            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    )
endif()

