#include <stdexcept>
#include <iostream>

#include <tclap/CmdLine.h>

#include <andres/functional.hxx>
#include <andres/lp/gurobi.hxx>
#include <andres/ilp/gurobi.hxx>
#include <andres/graph/multicut-lifted/ilp.hxx>
#include <andres/graph/multicut-lifted/lp.hxx>
#include <andres/graph/multicut-lifted/greedy-additive.hxx>
#include <andres/graph/multicut-lifted/kernighan-lin.hxx>

#include "Timer.hpp"
#include "utils.hxx"



enum class Method {
    Zeros,
    Ones,
    GAEC,
    Kernighan_Lin,
    LP,
    ILP
};

enum class Initialization {
    Zeros,
    Ones,
    Input_Labeling,
    GAEC
};

struct Parameters {
    std::string inputHDF5FileName;
    std::string outputHDF5FileName;
    std::string labelingHDF5FileName;
    Method optimizationMethod { Method::Kernighan_Lin };
    Initialization initialization { Initialization::Zeros };
    bool probabilistic { true };
};

inline void
parseCommandLine(
    int argc,
    char** argv,
    Parameters& parameters
)
try
{
    TCLAP::CmdLine tclap("solve-lifted-multicut-problem-grid-graph", ' ', "1.0");
    TCLAP::ValueArg<std::string> argInputHDF5FileName("i", "input-hdf5-file", "File to load multicut problem from", true, parameters.inputHDF5FileName, "INPUT_HDF5_FILE", tclap);
    TCLAP::ValueArg<std::string> argOutputHDF5FileName("o", "output-hdf-file", "hdf file (output)", false, parameters.outputHDF5FileName, "OUTPUT_HDF5_FILE", tclap);
    TCLAP::ValueArg<std::string> argLabelingHDF5FileName("l", "labeling-hdf-file", "hdf file specifying initial node labelings (input)", false, parameters.labelingHDF5FileName, "LABELING_HDF5_FILE", tclap);
    TCLAP::ValueArg<std::string> argOptimizationMethod("m", "optimization-method", "optimization method to use {zeros, ones, ILP, LP, GAEC, KL}", false, "KL", "OPTIMIZATION_METHOD", tclap);
    TCLAP::ValueArg<std::string> argInitializationMethod("I", "initialization-method", "initialization method to use {zeros, ones, GAEC}", false, "zeros", "INITIALIZATION_METHOD", tclap);
    TCLAP::SwitchArg argNonProbabilistic("p", "non-probabilistic", "Assume inputs are not probabilities. By default, all inputs are assumed to be Logistic Probabilities. (Default: disabled).",tclap);

    tclap.parse(argc, argv);

    parameters.inputHDF5FileName = argInputHDF5FileName.getValue();
    parameters.outputHDF5FileName = argOutputHDF5FileName.getValue();
    parameters.probabilistic = !argNonProbabilistic.getValue();

    if (!argOptimizationMethod.isSet())
        throw std::runtime_error("No optimization method specified");
    
    if (argOptimizationMethod.getValue() == "GAEC")
        parameters.optimizationMethod = Method::GAEC;
    else if (argOptimizationMethod.getValue() == "zeros")
        parameters.optimizationMethod = Method::Zeros;
    else if (argOptimizationMethod.getValue() == "ones")
        parameters.optimizationMethod = Method::Ones;
    else if (argOptimizationMethod.getValue() == "LP")
        parameters.optimizationMethod = Method::LP;
    else if (argOptimizationMethod.getValue() == "ILP")
        parameters.optimizationMethod = Method::ILP;
    else if(argOptimizationMethod.getValue() == "KL")
    {
        parameters.optimizationMethod = Method::Kernighan_Lin;

        if (!argInitializationMethod.isSet() && !argLabelingHDF5FileName.isSet())
            throw std::runtime_error("Either initialization method (zeros, ones) or initial labeling must be specified for Kernighan-Lin.");

        if(argLabelingHDF5FileName.isSet())
        {
            if(argInitializationMethod.isSet())
                throw std::runtime_error("Either initialization method or initial labeling must be specified.");
            
            parameters.labelingHDF5FileName = argLabelingHDF5FileName.getValue();
            parameters.initialization = Initialization::Input_Labeling;
        }
        else if (argInitializationMethod.isSet())
        {
            if(argInitializationMethod.getValue() == "ones")
                parameters.initialization = Initialization::Ones;
            else if(argInitializationMethod.getValue() == "zeros")
                parameters.initialization = Initialization::Zeros;
            else if(argInitializationMethod.getValue() == "GAEC")
                parameters.initialization = Initialization::GAEC;
            else
                throw std::runtime_error("Invalid initialization method specified");
        }
    }
    else
        throw std::runtime_error("Invalid optimization method specified");
}
catch(TCLAP::ArgException& e)
{
    throw std::runtime_error(e.error());
}

template<typename GraphType>
void solveLiftedMulticutProblem(
    const Parameters& parameters,
    std::ostream& stream = std::cerr
)
{
    GraphType original_graph;
    andres::graph::Graph<> lifted_graph;
    std::vector<double> edge_values;

    // Load Lifted Multicut Problem
    {
        auto fileHandle = andres::graph::hdf5::openFile(parameters.inputHDF5FileName);

        andres::graph::hdf5::load(fileHandle, "graph", original_graph);        
        andres::graph::hdf5::load(fileHandle, "graph-lifted", lifted_graph);

        std::vector<size_t> shape;
        andres::graph::hdf5::load(fileHandle, "edge-cut-probabilities", shape, edge_values);
        andres::graph::hdf5::closeFile(fileHandle);
    }

    if (parameters.probabilistic)
        std::transform(
            edge_values.begin(),
            edge_values.end(),
            edge_values.begin(),
            andres::NegativeLogProbabilityRatio<double,double>()
            );

    // Solve Lifted Multicut problem
    std::vector<char> edge_labels(lifted_graph.numberOfEdges());
    
    Timer t;
    t.start();

    if (parameters.initialization == Initialization::Zeros)
        std::fill(edge_labels.begin(), edge_labels.end(), 0);
    else if (parameters.initialization == Initialization::Ones)
        std::fill(edge_labels.begin(), edge_labels.end(), 1);
    else if (parameters.initialization == Initialization::Input_Labeling)
    {
        auto fileHandle = andres::graph::hdf5::openFile(parameters.labelingHDF5FileName);

        std::vector<size_t> shape;
        std::vector<size_t> vertex_labels;
        andres::graph::hdf5::load(fileHandle, "labels", shape, vertex_labels);

        andres::graph::hdf5::closeFile(fileHandle);

        assert(shape.size() == 1);
        assert(shape[0] == lifted_graph.numberOfVertices());

        edge_labels.resize(lifted_graph.numberOfEdges());

        vertexToEdgeLabels(original_graph, lifted_graph, vertex_labels, edge_labels);
    }
    else if (parameters.initialization == Initialization::GAEC)
        andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(original_graph, lifted_graph, edge_values, edge_labels);

    if (parameters.optimizationMethod == Method::Ones)
        std::fill(edge_labels.begin(), edge_labels.end(), 1);
    else if (parameters.optimizationMethod == Method::Zeros)
        std::fill(edge_labels.begin(), edge_labels.end(), 0);
    else if (parameters.optimizationMethod == Method::ILP)
        andres::graph::multicut_lifted::ilp<andres::ilp::Gurobi<>>(original_graph, lifted_graph, edge_values, edge_labels, edge_labels);
    else if (parameters.optimizationMethod == Method::GAEC)
        andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(original_graph, lifted_graph, edge_values, edge_labels);
    else if (parameters.optimizationMethod == Method::Kernighan_Lin)
        andres::graph::multicut_lifted::kernighanLin(original_graph, lifted_graph, edge_values, edge_labels, edge_labels);
    else if (parameters.optimizationMethod == Method::LP)
    {
        auto values = andres::graph::multicut_lifted::lp<andres::relax::Gurobi<>>(original_graph, lifted_graph, edge_values);

        t.stop();

        auto energy_value = inner_product(edge_values.begin(), edge_values.end(), values.begin(), .0);

        if (!parameters.outputHDF5FileName.empty())
        {
            auto file = andres::graph::hdf5::createFile(parameters.outputHDF5FileName);
            
            andres::graph::hdf5::save(file, "graph", original_graph);
            andres::graph::hdf5::save(file, "energy-value", energy_value);
            andres::graph::hdf5::save(file, "running-time", t.get_elapsed_seconds());
            andres::graph::hdf5::save(file, "labels", { values.size() }, values.data()); // we save directly edge values as given by the LP solution, since the latter is not guaranteed to be integer

            andres::graph::hdf5::closeFile(file);
        }

        std::cout << "Number of clusters: N/A\n";
        std::cout << "Energy value: " << energy_value << std::endl;
        std::cout << "Running time: " << t.to_string() << std::endl;

        return;
    }
    else
        throw std::runtime_error("Unsupported algorithm");
    
    t.stop();
    
    stream << "saving decomposition into file: " << parameters.outputHDF5FileName << std::endl;
    {
        auto file = andres::graph::hdf5::createFile(parameters.outputHDF5FileName);
        
        andres::graph::hdf5::save(file, "graph", original_graph);

        std::vector<size_t> vertex_labels(lifted_graph.numberOfVertices());
        edgeToVertexLabels(lifted_graph, edge_labels, vertex_labels);

        andres::graph::hdf5::save(file, "labels", { vertex_labels.size() }, vertex_labels.data());

        auto energy_value = inner_product(edge_values.begin(), edge_values.end(), edge_labels.begin(), .0);

        andres::graph::hdf5::save(file, "energy-value", energy_value);
        andres::graph::hdf5::save(file, "running-time", t.get_elapsed_seconds());

        std::vector<char> true_edge_labels(lifted_graph.numberOfEdges());
        vertexToEdgeLabels(original_graph, lifted_graph, vertex_labels, true_edge_labels);

        auto true_energy_value = inner_product(edge_values.begin(), edge_values.end(), true_edge_labels.begin(), .0);        

        andres::graph::hdf5::save(file, "true-energy-value", true_energy_value);

        andres::graph::hdf5::closeFile(file);

        std::cout << "Number of clusters: " << *max_element(vertex_labels.begin(), vertex_labels.end()) + 1 << std::endl;
        std::cout << "Energy value: " << energy_value << std::endl;
        std::cout << "Running time: " << t.to_string() << std::endl;
        std::cout << "True energy value: " << true_energy_value << std::endl;
    }
}