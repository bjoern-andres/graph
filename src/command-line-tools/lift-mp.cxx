#include <stdexcept>
#include <iostream>
#include <sstream>
#include <numeric>

#include <tclap/CmdLine.h>

#include <andres/functional.hxx>
#include <andres/graph/graph.hxx>
#include <andres/graph/hdf5/graph.hxx>
#include <andres/graph/lifting.hxx>

#include "Timer.hpp"



using namespace std;
using namespace andres;
using namespace andres::graph;

typedef double value_type;
typedef size_t size_type;
typedef vector<value_type> ValueMap;

struct Parameters {
    string inputHDF5FileName_;
    string outputHDF5FileName_;
    size_type upperDistanceBound_ { 1 };
    size_type lowerDistanceBound_ { 0 };
    bool probabilistic_ { true };
};

inline void
parseCommandLine(
    int argc,
    char** argv,
    Parameters& parameters
)
try
{
    TCLAP::CmdLine tclap("lift-multicut-problem", ' ', "1.0");
    TCLAP::ValueArg<string> argInputHDF5FileName("i", "input-hdf5-file", "File to load multicut problem from", true, parameters.inputHDF5FileName_, "INPUT_HDF5_FILE",tclap);
    TCLAP::ValueArg<string> argOutputHDF5FileName("o", "output-hdf5-file", "hdf file (output)", true, parameters.outputHDF5FileName_, "OUTPUT_HDF5_FILE",tclap);
    TCLAP::ValueArg<size_t> argUpperDistanceBound("u", "distance-upper-bound", "upper distance bound (inclusive)", false, parameters.upperDistanceBound_, "DIST_UPPER",tclap);
    TCLAP::ValueArg<size_t> argLowerDistanceBound("l", "distance-lower-bound", "lower distance bound (non-inclusive)", false, parameters.lowerDistanceBound_, "DIST_LOWER",tclap);
    TCLAP::SwitchArg argNonProbabilistic("p", "non-probabilistic", "Assume inputs are not probabilities. By default, all inputs are assumed to be Logistic Probabilities. (Default: disabled).",tclap);

    tclap.parse(argc, argv);

    parameters.inputHDF5FileName_ = argInputHDF5FileName.getValue();
    parameters.outputHDF5FileName_ = argOutputHDF5FileName.getValue();
    parameters.upperDistanceBound_ = argUpperDistanceBound.getValue();
    parameters.probabilistic_ = !argNonProbabilistic.getValue();
    parameters.lowerDistanceBound_ = argLowerDistanceBound.getValue();
    parameters.upperDistanceBound_ = argUpperDistanceBound.getValue();
}
catch(TCLAP::ArgException& e) {
    throw runtime_error(e.error());
}

void liftMulticutProblem(
    const Parameters& parameters,
    ostream& stream = cerr
) { 
    {
        stream << "Liftin using neighbourhood defined by bounds (" << parameters.lowerDistanceBound_ 
            << ",";

        if (parameters.upperDistanceBound_ == numeric_limits<size_type>::max())
            stream << "MAX";
        else
            stream << parameters.upperDistanceBound_;

        stream << ")" << endl;
    }
    stream << "assuming Graph input." << endl;
    
    Graph<> graph;
    ValueMap edgeCutProbabilities;

    stream << "loading multicut problem from file: " << parameters.inputHDF5FileName_ <<  endl;
    {
        auto fileHandle = hdf5::openFile(parameters.inputHDF5FileName_);
        
        hdf5::load(fileHandle, "graph", graph);

        vector<size_t> shape;
        hdf5::load(fileHandle, "edge-cut-probabilities", shape, edgeCutProbabilities);

        hdf5::closeFile(fileHandle);

        assert(shape.size() == 1);
        assert(shape[0] == graph.numberOfEdges());
    }

    Timer t;
    t.start();
    
    Graph<> graphLifted;
    lift(graph, graphLifted, parameters.upperDistanceBound_, parameters.lowerDistanceBound_);
    
    ValueMap edgeSplitProbabilitiesLifted(graphLifted.numberOfEdges());
    if (parameters.probabilistic_)
        transform(
            edgeCutProbabilities.begin(),
            edgeCutProbabilities.end(),
            edgeCutProbabilities.begin(),
            ProbabilityToNegativeLogInverseProbability<value_type,value_type>()
        );
    
    liftEdgeValues(
        graph,
        graphLifted,
        edgeCutProbabilities.begin(),
        edgeSplitProbabilitiesLifted.begin()
    );
    
    if (parameters.probabilistic_)
        transform(
            edgeSplitProbabilitiesLifted.begin(),
            edgeSplitProbabilitiesLifted.end(),
            edgeSplitProbabilitiesLifted.begin(),
            NegativeLogProbabilityToInverseProbability<value_type,value_type>()
        );

    t.stop();
    
    stream << "saving Lifted Multicut as HDF5 file: " << parameters.outputHDF5FileName_ << endl;
    auto file = hdf5::createFile(parameters.outputHDF5FileName_);
    hdf5::save(file, "graph", graph);
    hdf5::save(file, "graph-lifted", graphLifted);
    hdf5::save(file, "lifting-time", t.get_elapsed_seconds());
    hdf5::save(file, "edge-cut-probabilities", { graphLifted.numberOfEdges() }, edgeSplitProbabilitiesLifted.data());
    hdf5::closeFile(file);
}

int main(int argc, char** argv)
try
{
    Parameters parameters;
    
    parseCommandLine(argc, argv, parameters);
    liftMulticutProblem(parameters);

    return 0;
}
catch(const runtime_error& error)
{
    cerr << "error creating multicut problem: " << error.what() << endl;
    return 1;
}