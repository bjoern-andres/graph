#pragma once

#include <stdexcept>
#include <iostream>
#include <numeric>

#include <tclap/CmdLine.h>

#include <andres/ilp/gurobi.hxx>
#include <andres/lp/gurobi.hxx>

#include <andres/graph/multicut/ilp.hxx>
#include <andres/graph/multicut/lp.hxx>
#include <andres/graph/multicut/ilp-callback.hxx>
#include <andres/graph/multicut/greedy-additive.hxx>
#include <andres/graph/multicut/greedy-fixation.hxx>
#include <andres/graph/multicut/kernighan-lin.hxx>

#include "Timer.hpp"
#include "utils.hxx"



enum class Method {
    Zeros,
    Ones,
    GAEC,
    GF,
    KL,
    ILP,
    ILPC,
    LP
};

enum class Initialization {
    Zeros,
    Ones,
    Input_Labeling,
    GAEC,
    GF
};

struct Parameters {
    string inputHDF5FileName;
    string outputHDF5FileName;
    string labelingHDF5FileName;
    Method optimizationMethod_;
    Initialization initialization;
};

inline void
parseCommandLine(
    int argc,
    char** argv,
    Parameters& parameters
)
{
    try
    {
        TCLAP::CmdLine tclap("solve-multicut-problem", ' ', "1.0");
        TCLAP::ValueArg<string> argInputHDF5FileName("i", "hdf5-input", "File to load multicut problem from", true, parameters.inputHDF5FileName, "INPUT_HDF5_FILE", tclap);
        TCLAP::ValueArg<string> argOutputHDF5FileName("o", "output-hdf-file", "hdf file (output)", true, parameters.outputHDF5FileName, "OUTPUT_HDF5_FILE", tclap);
        TCLAP::ValueArg<string> argLabelingHDF5FileName("l", "labeling-hdf-file", "hdf file specifying initial node labelings (input)", false, parameters.labelingHDF5FileName, "LABELING_HDF5_FILE", tclap);

        TCLAP::ValueArg<string> argOptimizationMethod("m", "optimization-method", "optimization method to use {LP, ILP, ILPC, GAEC, GF, KL, zeros, ones}", false, "KL", "OPTIMIZATION_METHOD", tclap);
        TCLAP::ValueArg<string> argInitializationMethod("I", "initialization-method", "initialization method to use {zeros, ones, GAEC, GF}", false, "zeros", "INITIALIZATION_METHOD", tclap);
        
        tclap.parse(argc, argv);

        parameters.inputHDF5FileName = argInputHDF5FileName.getValue();
        parameters.outputHDF5FileName = argOutputHDF5FileName.getValue();
        parameters.labelingHDF5FileName = argLabelingHDF5FileName.getValue();

        if (!argOptimizationMethod.isSet())
            throw std::runtime_error("No optimization method specified");

        if (!parameters.labelingHDF5FileName.empty())
            parameters.initialization = Initialization::Input_Labeling;
        else if(argInitializationMethod.getValue() == "ones")
            parameters.initialization = Initialization::Ones;
        else if(argInitializationMethod.getValue() == "zeros")
            parameters.initialization = Initialization::Zeros;
        else if(argInitializationMethod.getValue() == "GAEC")
            parameters.initialization = Initialization::GAEC;
        else if(argInitializationMethod.getValue() == "GF")
            parameters.initialization = Initialization::GF;
        else
            throw std::runtime_error("Invalid initialization method specified");
        
        if (argOptimizationMethod.getValue() == "GAEC")
            parameters.optimizationMethod_ = Method::GAEC;
        else if(argOptimizationMethod.getValue() == "GF")
            parameters.optimizationMethod_ = Method::GF;
        else if(argOptimizationMethod.getValue() == "KL")
            parameters.optimizationMethod_ = Method::KL;
        else if (argOptimizationMethod.getValue() == "ILP")
            parameters.optimizationMethod_ = Method::ILP;
        else if (argOptimizationMethod.getValue() == "ILPC")
            parameters.optimizationMethod_ = Method::ILPC;
        else if (argOptimizationMethod.getValue() == "LP")
            parameters.optimizationMethod_ = Method::LP;
        else if(argInitializationMethod.getValue() == "ones")
            parameters.optimizationMethod_ = Method::Ones;
        else if(argInitializationMethod.getValue() == "zeros")
            parameters.optimizationMethod_ = Method::Zeros;
        else
            throw std::runtime_error("Invalid optimization method specified");
    }
    catch(TCLAP::ArgException& e)
    {
        throw std::runtime_error(e.error());
    }
}

template<typename GraphType>
void solveMulticutProblem(
    const Parameters& parameters
)
{
    GraphType graph;
    std::vector<double> edge_values;

    {
        auto fileHandle = andres::graph::hdf5::openFile(parameters.inputHDF5FileName);

        andres::graph::hdf5::load(fileHandle, "graph", graph);

        std::vector<size_t> shape;
        andres::graph::hdf5::load(fileHandle, "edge-values", shape, edge_values);

        andres::graph::hdf5::closeFile(fileHandle);
    }
    std::cout << "Number of nodes: " << graph.numberOfVertices() << std::endl;
    std::cout << "Number of edges: " << graph.numberOfEdges() << std::endl;

    // Solve Multicut problem
    std::vector<char> edge_labels(graph.numberOfEdges());

    Timer t;
    t.start();

    if (parameters.initialization == Initialization::Zeros)
        fill(edge_labels.begin(), edge_labels.end(), 0);
    else if (parameters.initialization == Initialization::Ones)
        fill(edge_labels.begin(), edge_labels.end(), 1);
    else if (parameters.initialization == Initialization::GAEC)
        andres::graph::multicut::greedyAdditiveEdgeContraction(graph, edge_values, edge_labels);
    else if (parameters.initialization == Initialization::GF)
        andres::graph::multicut::greedyFixation(graph, edge_values, edge_labels);
    else if (parameters.initialization == Initialization::Input_Labeling)
    {
        auto fileHandle = andres::graph::hdf5::openFile(parameters.labelingHDF5FileName);

        std::vector<size_t> shape;
        std::vector<size_t> vertex_labels;
        andres::graph::hdf5::load(fileHandle, "labels", shape, vertex_labels);

        andres::graph::hdf5::closeFile(fileHandle);

        assert(shape.size() == 1);
        assert(shape[0] == graph.numberOfVertices());

        edge_labels.resize(graph.numberOfEdges());
        vertexToEdgeLabels(graph, vertex_labels, edge_labels);
    }

    if (parameters.optimizationMethod_ == Method::Zeros)
        fill(edge_labels.begin(), edge_labels.end(), 0);
    else if (parameters.optimizationMethod_ == Method::Ones)
        fill(edge_labels.begin(), edge_labels.end(), 1);
    else if (parameters.optimizationMethod_ == Method::GAEC)
        andres::graph::multicut::greedyAdditiveEdgeContraction(graph, edge_values,  edge_labels);
    else if (parameters.optimizationMethod_ == Method::GF)
        andres::graph::multicut::greedyFixation(graph, edge_values,  edge_labels);
    else if (parameters.optimizationMethod_ == Method::KL)
        andres::graph::multicut::kernighanLin(graph, edge_values, edge_labels, edge_labels);
    else if (parameters.optimizationMethod_ == Method::ILP)
        andres::graph::multicut::ilp<andres::ilp::Gurobi>(graph, edge_values, edge_labels, edge_labels);
    else if (parameters.optimizationMethod_ == Method::ILPC)
        andres::graph::multicut::ilp<andres::ilp::Gurobi>(graph, edge_values, edge_labels, edge_labels);
    else if (parameters.optimizationMethod_ == Method::LP)
    {
        auto values = andres::graph::multicut::lp<andres::lp::Gurobi>(graph, edge_values);

        t.stop();

        auto energy_value = inner_product(edge_values.begin(), edge_values.end(), values.begin(), .0);

        if (!parameters.outputHDF5FileName.empty())
        {
            auto file = andres::graph::hdf5::createFile(parameters.outputHDF5FileName);
            
            andres::graph::hdf5::save(file, "graph", graph);
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
    
    auto energy_value = inner_product(edge_values.begin(), edge_values.end(), edge_labels.begin(), .0);

    std::vector<size_t> vertex_labels(graph.numberOfVertices());
    edgeToVertexLabels(graph, edge_labels, vertex_labels);

    auto file = andres::graph::hdf5::createFile(parameters.outputHDF5FileName);

    andres::graph::hdf5::save(file, "graph", graph);
    andres::graph::hdf5::save(file, "labels", { vertex_labels.size() }, vertex_labels.data());
    andres::graph::hdf5::save(file, "energy-value", energy_value);
    andres::graph::hdf5::save(file, "running-time", t.get_elapsed_seconds());
    andres::graph::hdf5::closeFile(file);
    
    std::cout << "Number of clusters: " << *max_element(vertex_labels.begin(), vertex_labels.end()) + 1 << std::endl;
    std::cout << "Energy value: " << energy_value << std::endl;
    std::cout << "Running time: " << t.to_string() << std::endl;
}