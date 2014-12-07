
#pragma once
#ifndef ANDRES_GRAPH_COMPLETE_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_COMPLETE_GRAPH_HDF5_HXX

#include <array>

#include "hdf5.hxx"
#include "../complete-graph.hxx"


namespace andres{
namespace graph{
namespace hdf5{


struct DatasetNamesComplete : public DatasetNames {
    const char vertices[9] = "vertices";
};
static const DatasetNamesComplete datasetNamesComplete;

template<>
template<class VISITOR>
struct graph_traits<CompleteGraph<VISITOR> > {
    static const int ID = 10002;
};

template <class VISITOR>
void save(const hid_t, const std::string&, const CompleteGraph<VISITOR>& graph);

template <class VISITOR>
void load(const hid_t, const std::string&, CompleteGraph<VISITOR>& graph);

// Implementation

template <class VISITOR>
void
save(
    const hid_t fileHandle,
    const std::string& datasetName,
    const CompleteGraph<VISITOR>& graph
) {
    typedef CompleteGraph<VISITOR> Graph;
    hdf5::HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(fileHandle, datasetName,true);
    std::string sError;
    try{
        save(groupHandle, datasetNamesComplete.vertices, graph.numberOfVertices());
        int ID = graph_traits<Graph>::ID;
        save(groupHandle, datasetNamesComplete.graphTypeId, ID);
    } catch (std::exception& e) {
        sError = e.what();
    }
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "CompleteGraph: Saving to dataset '"+datasetName+"' failed: "+sError
        );
    }
}

template <class VISITOR>
void
load(
    const hid_t fileHandle,
    const std::string& datasetName,
    CompleteGraph<VISITOR>& graph
) {
    typedef CompleteGraph<VISITOR> CompleteGraph;
    hdf5::HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(fileHandle, datasetName);
    std::string sError;
    try{
        {
            int ID;
            load(groupHandle, datasetNamesComplete.graphTypeId, ID);
            if(ID!=graph_traits<CompleteGraph>::ID) {
                sError = "Stored graph type is not a CompleteGraph.";
                goto cleanup;
            }
        }
        std::size_t numberOfVertices;
        load(groupHandle, datasetNamesComplete.vertices, numberOfVertices);
        graph.assign(numberOfVertices);
        }catch(std::exception& e) {
            sError = e.what();
    }
cleanup:
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "CompleteGraph: Loading dataset '"+datasetName+"' failed: "+sError
        );
    }
}


} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_COMPLETE_GRAPH_HDF5_HXX