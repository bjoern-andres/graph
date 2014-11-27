cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(graph)
include_directories(include)
file(GLOB headers
    include/andres/*.hxx
    include/andres/graph/*.hxx
    include/andres/graph/detail/*.hxx
    include/andres/ilp/*.hxx
)
enable_testing()

##############################################################################
# Doxygen
##############################################################################
find_package(Doxygen)

##############################################################################
# OpenMP
##############################################################################
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

##############################################################################
# Gurobi
##############################################################################
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIR})
endif()

##############################################################################
# Cplex
##############################################################################
find_package(CPLEX)
if(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    include_directories(${CPLEX_CONCERT_INCLUDE_DIR})
    add_definitions(-DILOUSESTL)
    add_definitions(-DIL_STD)
endif()

##############################################################################
# C++11 support
##############################################################################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Some functionality will not be available.")
endif()

##############################################################################
# targets: graph
##############################################################################
add_executable(test-graph-graph src/andres/graph/unit-test/graph.cxx ${headers})
add_test(test-graph-graph test-graph-graph)

add_executable(test-graph-digraph src/andres/graph/unit-test/digraph.cxx ${headers})
add_test(test-graph-digraph test-graph-digraph)

add_executable(test-graph-complete-graph src/andres/graph/unit-test/graph-complete.cxx ${headers})
add_test(test-graph-complete-graph test-graph-complete-graph)

if(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)
    add_executable(test-graph-grid-graph src/andres/graph/unit-test/graph-grid.cxx ${headers})
    add_test(test-graph-grid-graph test-graph-grid-graph)
endif(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)

add_executable(test-graph-shortest-paths src/andres/graph/unit-test/shortest-paths.cxx ${headers})
add_test(test-graph-shortest-paths test-graph-shortest-paths)

add_executable(test-graph-paths src/andres/graph/unit-test/paths.cxx ${headers})
add_test(test-graph-paths test-graph-paths)

add_executable(test-graph-components src/andres/graph/unit-test/components.cxx ${headers})
add_test(test-graph-components test-graph-components)

add_executable(test-graph-max-flow src/andres/graph/unit-test/max-flow.cxx ${headers})
add_test(test-graph-max-flow test-graph-max-flow)

if(CPLEX_FOUND)
    add_executable(test-graph-multicut-cplex src/andres/graph/unit-test/multicut-cplex.cxx ${headers})
    target_link_libraries(test-graph-multicut-cplex ${CPLEX_LIBRARIES})
    add_test(test-graph-multicut-cplex test-graph-multicut-cplex)
endif()

if(GUROBI_FOUND)
    add_executable(test-graph-multicut-gurobi src/andres/graph/unit-test/multicut-gurobi.cxx ${headers})
    target_link_libraries(test-graph-multicut-gurobi ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-gurobi test-graph-multicut-gurobi)
endif()

if(DOXYGEN_FOUND)
    configure_file("${graph_SOURCE_DIR}/doxygen/doxyfile-graph.in" "${graph_BINARY_DIR}/doxyfile-graph" @ONLY IMMEDIATE)
    add_custom_target(doc-graph ALL COMMAND ${DOXYGEN} "${graph_BINARY_DIR}/doxyfile-graph")
endif()
