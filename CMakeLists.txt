cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(graph)
include_directories(include)
file(GLOB headers include/andres/*.hxx include/andres/graph/*.hxx include/andres/ilp/*.hxx)
enable_testing()

##############################################################################
# Doxygen
##############################################################################
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    message(STATUS "Doxygen found")
else()
    message("doxygen not found")
endif()

##############################################################################
# OpenMP
##############################################################################
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found")
else()
    message("OpenMP not found")
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
    include_directories(${CONCERT_INCLUDE_DIR})
    add_definitions(-DILOUSESTL)
    add_definitions(-DIL_STD)
endif()

##############################################################################
# targets: graph
##############################################################################
add_executable(test-graph-graph src/unittest/graph/graph.cxx ${headers})
add_test(test-graph-graph test-graph-graph)

add_executable(test-graph-digraph src/unittest/graph/digraph.cxx ${headers})
add_test(test-graph-digraph test-graph-digraph)

add_executable(test-graph-shortest-paths src/unittest/graph/shortest-paths.cxx ${headers})
add_test(test-graph-shortest-paths test-graph-shortest-paths)

add_executable(test-graph-paths src/unittest/graph/paths.cxx ${headers})
add_test(test-graph-paths test-graph-paths)

add_executable(test-graph-components src/unittest/graph/components.cxx ${headers})
add_test(test-graph-components test-graph-components)

add_executable(test-max-flow src/unittest/graph/max-flow.cxx ${headers})
add_test(test-max-flow test-max-flow)

if(CPLEX_FOUND)
    add_executable(test-graph-multicut-cplex src/unittest/graph/multicut-cplex.cxx ${headers})
    target_link_libraries(test-graph-multicut-cplex ${CMAKE_THREAD_LIBS_INIT} ${ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY} ${CONCERT_LIBRARY})
    add_test(test-graph-multicut-cplex test-graph-multicut-cplex)
endif()

if(GUROBI_FOUND)
    add_executable(test-graph-multicut-gurobi src/unittest/graph/multicut-gurobi.cxx ${headers})
    target_link_libraries(test-graph-multicut-gurobi ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-gurobi test-graph-multicut-gurobi)
endif()

if(DOXYGEN_FOUND)
    configure_file("${graph_SOURCE_DIR}/doxygen/doxyfile-graph.in" "${graph_BINARY_DIR}/doxyfile-graph" @ONLY IMMEDIATE)
    add_custom_target(doc-graph ALL COMMAND ${DOXYGEN} "${graph_BINARY_DIR}/doxyfile-graph")
endif()
