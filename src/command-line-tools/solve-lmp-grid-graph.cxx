#include <andres/graph/graph.hxx>
#include <andres/graph/hdf5/graph.hxx>
#include <andres/graph/grid-graph.hxx>
#include <andres/graph/hdf5/grid-graph.hxx>

#include "solve-lmp.hxx"

int main(int argc, char** argv)
try 
{
    Parameters parameters;

    parseCommandLine(argc, argv, parameters);
    solveLiftedMulticutProblem<andres::graph::GridGraph<>>(parameters);

    return 0;
}
catch (const std::runtime_error& error)
{
    std::cerr << "error creating multicut problem: " << error.what() << std::endl;
    return 1;
}