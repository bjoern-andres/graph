cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(graph)
enable_testing()

##############################################################################
# Doxygen
##############################################################################
find_package(Doxygen)

##############################################################################
# OpenMP
##############################################################################
#find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

##############################################################################
# HDF5
##############################################################################
find_package(HDF5 1.8.0)
include_directories(${HDF5_INCLUDE_DIR})

##############################################################################
# Gurobi
##############################################################################
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIR})
    add_definitions("-DWITH_GUROBI")
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
add_executable(test-graph-graph src/andres/graph/unit-test/graph.cxx)
target_include_directories(test-graph-graph PUBLIC include)
add_test(test-graph-graph test-graph-graph)

add_executable(test-graph-digraph src/andres/graph/unit-test/digraph.cxx)
target_include_directories(test-graph-digraph PUBLIC include)
add_test(test-graph-digraph test-graph-digraph)

add_executable(test-graph-shortest-paths src/andres/graph/unit-test/shortest-paths.cxx)
target_include_directories(test-graph-shortest-paths PUBLIC include)
add_test(test-graph-shortest-paths test-graph-shortest-paths)

add_executable(test-graph-paths src/andres/graph/unit-test/paths.cxx)
target_include_directories(test-graph-paths PUBLIC include)
add_test(test-graph-paths test-graph-paths)

add_executable(test-graph-components src/andres/graph/unit-test/components.cxx)
target_include_directories(test-graph-components PUBLIC include)
add_test(test-graph-components test-graph-components)

add_executable(test-graph-max-flow src/andres/graph/unit-test/max-flow.cxx)
target_include_directories(test-graph-max-flow PUBLIC include)
add_test(test-graph-max-flow test-graph-max-flow)

add_executable(analysis-erase-vertex src/andres/graph/analysis/erase-vertex.cxx)
target_include_directories(analysis-erase-vertex PUBLIC include)

if(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)
    add_executable(test-graph-complete-graph src/andres/graph/unit-test/graph-complete.cxx)
    target_include_directories(test-graph-complete-graph PUBLIC include)
    add_test(test-graph-complete-graph test-graph-complete-graph)

    add_executable(test-graph-grid-graph src/andres/graph/unit-test/graph-grid.cxx)
    target_include_directories(test-graph-grid-graph PUBLIC include)
    add_test(test-graph-grid-graph test-graph-grid-graph)

    add_executable(test-graph-dfs src/andres/graph/unit-test/dfs.cxx)
    target_include_directories(test-graph-dfs PUBLIC include)
    add_test(test-graph-dfs test-graph-dfs)

    add_executable(test-graph-bfs src/andres/graph/unit-test/bfs.cxx)
    target_include_directories(test-graph-bfs PUBLIC include)
    add_test(test-graph-bfs test-graph-bfs)

    add_executable(test-graph-cut-vertices src/andres/graph/unit-test/cut-vertices.cxx)
    target_include_directories(test-graph-cut-vertices PUBLIC include)
    add_test(test-graph-cut-vertices test-graph-cut-vertices)

    add_executable(test-graph-bridges src/andres/graph/unit-test/bridges.cxx)
    target_include_directories(test-graph-bridges PUBLIC include)
    add_test(test-graph-bridges test-graph-bridges)

    add_executable(test-graph-triangles src/andres/graph/unit-test/triangles.cxx)
    target_include_directories(test-graph-triangles PUBLIC include)
    add_test(test-graph-triangles test-graph-triangles)

    add_executable(test-graph-minimum-spanning-tree src/andres/graph/unit-test/minimum-spanning-tree.cxx)
    target_include_directories(test-graph-minimum-spanning-tree PUBLIC include)
    add_test(test-graph-minimum-spanning-tree test-graph-minimum-spanning-tree)

    add_executable(test-graph-lifting src/andres/graph/unit-test/lifting.cxx)
    target_include_directories(test-graph-lifting PUBLIC include)
    add_test(test-graph-lifting test-graph-lifting)

    add_executable(test-graph-multicut-lifted-kl src/andres/graph/unit-test/multicut-lifted/kernighan-lin.cxx)
    target_include_directories(test-graph-multicut-lifted-kl PUBLIC include)
    add_test(test-graph-multicut-lifted-kl test-graph-multicut-lifted-kl)

    add_executable(test-graph-multicut-lifted-greedy-additive src/andres/graph/unit-test/multicut-lifted/greedy-additive.cxx)
    target_include_directories(test-graph-multicut-lifted-greedy-additive PUBLIC include)
    add_test(test-graph-multicut-lifted-greedy-additive test-graph-multicut-lifted-greedy-additive)

    add_executable(test-graph-multicut-kernighan-lin src/andres/graph/unit-test/multicut/kernighan-lin.cxx)
    target_include_directories(test-graph-multicut-kernighan-lin PUBLIC include)
    add_test(test-graph-multicut-kernighan-lin test-graph-multicut-kernighan-lin)

    add_executable(test-graph-multicut-greedy-additive src/andres/graph/unit-test/multicut/greedy-additive.cxx)
    target_include_directories(test-graph-multicut-greedy-additive PUBLIC include)
    add_test(test-graph-multicut-greedy-additive test-graph-multicut-greedy-additive)

    add_executable(test-graph-bipartite-matching src/andres/graph/unit-test/bipartite-matching.cxx)
    target_include_directories(test-graph-bipartite-matching PUBLIC include)
    add_test(test-graph-bipartite-matching test-graph-bipartite-matching)
    
    add_executable(test-graph-multicut-mutex-watershed src/andres/graph/unit-test/multicut/mutex-watershed.cxx)
    target_include_directories(test-graph-multicut-mutex-watershed PUBLIC include)
    add_test(test-graph-multicut-mutex-watersehd test-graph-multicut-mutex-watershed)

    add_executable(test-graph-multicut-preprocessing src/andres/graph/unit-test/multicut/preprocessing.cxx)
    target_include_directories(test-graph-multicut-preprocessing PUBLIC include)
    add_test(test-graph-multicut-preprocessing test-graph-multicut-preprocessing)

    if(HDF5_FOUND)
        add_executable(test-hdf5 src/andres/graph/unit-test/hdf5.cxx)
        target_include_directories(test-hdf5 PUBLIC include ${HDF5_INCLUDE_DIRS})
        target_link_libraries(test-hdf5 ${HDF5_LIBRARIES})
        add_test(test-hdf5 test-hdf5)
    endif()
endif(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)

if(GUROBI_FOUND)
    add_executable(test-graph-multicut-ilp src/andres/graph/unit-test/multicut/ilp.cxx)
    target_include_directories(test-graph-multicut-ilp PUBLIC include)
    target_link_libraries(test-graph-multicut-ilp ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-ilp test-graph-multicut-ilp)

    add_executable(test-graph-multicut-ilp-callback src/andres/graph/unit-test/multicut/ilp-callback.cxx)
    target_include_directories(test-graph-multicut-ilp-callback PUBLIC include)
    target_link_libraries(test-graph-multicut-ilp-callback ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-ilp-callback test-graph-multicut-ilp-callback)

    add_executable(test-graph-multicut-lifted-ilp src/andres/graph/unit-test/multicut-lifted/ilp.cxx)
    target_include_directories(test-graph-multicut-lifted-ilp PUBLIC include)
    target_link_libraries(test-graph-multicut-lifted-ilp ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-lifted-ilp test-graph-multicut-lifted-ilp)

    add_executable(test-graph-multicut-lifted-ilp-callback src/andres/graph/unit-test/multicut-lifted/ilp-callback.cxx)
    target_include_directories(test-graph-multicut-lifted-ilp-callback PUBLIC include)
    target_link_libraries(test-graph-multicut-lifted-ilp-callback ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    add_test(test-graph-multicut-lifted-ilp-callback test-graph-multicut-lifted-ilp-callback)
endif()

if(DOXYGEN_FOUND)
    configure_file("${graph_SOURCE_DIR}/doxygen/doxyfile-graph.in" "${graph_BINARY_DIR}/doxyfile-graph" @ONLY IMMEDIATE)
    add_custom_target(doc-graph ALL COMMAND ${DOXYGEN} "${graph_BINARY_DIR}/doxyfile-graph")
endif()

if((COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11) AND HDF5_FOUND)
    add_executable(test-probabilistic-lifting src/command-line-tools/test-probabilistic-lifting.cxx)
    target_include_directories(test-probabilistic-lifting PUBLIC include src)
    add_test(test-probabilistic-lifting test-probabilistic-lifting)

    add_executable(lift-mp src/command-line-tools/lift-mp.cxx)
    target_include_directories(lift-mp PUBLIC include src)
    target_link_libraries(lift-mp ${HDF5_LIBRARIES})

    add_executable(lift-mp-grid-graph src/command-line-tools/lift-mp-grid-graph.cxx)
    target_include_directories(lift-mp-grid-graph PUBLIC include src)
    target_link_libraries(lift-mp-grid-graph ${HDF5_LIBRARIES})

    add_executable(solve-mp src/command-line-tools/solve-mp.cxx)
    target_include_directories(solve-mp PUBLIC include src)
    target_link_libraries(solve-mp ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})

    add_executable(solve-mp-complete-graph src/command-line-tools/solve-mp-complete-graph.cxx)
    target_include_directories(solve-mp-complete-graph PUBLIC include src)
    target_link_libraries(solve-mp-complete-graph ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})

    add_executable(solve-mp-grid-graph src/command-line-tools/solve-mp-grid-graph.cxx)
    target_include_directories(solve-mp-grid-graph PUBLIC include src)
    target_link_libraries(solve-mp-grid-graph ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})

    add_executable(solve-lmp src/command-line-tools/solve-lmp.cxx)
    target_include_directories(solve-lmp PUBLIC include src)
    target_link_libraries(solve-lmp ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})

    add_executable(solve-lmp-grid-graph src/command-line-tools/solve-lmp-grid-graph.cxx)
    target_include_directories(solve-lmp-grid-graph PUBLIC include src)
    target_link_libraries(solve-lmp-grid-graph ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})
endif()
