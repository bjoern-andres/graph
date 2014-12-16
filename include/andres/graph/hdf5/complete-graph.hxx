
#pragma once
#ifndef ANDRES_GRAPH_COMPLETE_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_COMPLETE_GRAPH_HDF5_HXX

#include <array>

#include "hdf5.hxx"
#include "../complete-graph.hxx"


namespace andres{
namespace graph{
namespace hdf5{


template<>
template<class VISITOR>
struct GraphTraitsHDF5<CompleteGraph<VISITOR> > {
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
    const hid_t parentHandle,
    const std::string& graphName,
    const CompleteGraph<VISITOR>& graph
) {
    typedef CompleteGraph<VISITOR> Graph;
    hdf5::HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName,true);
    std::string sError;
    try{
        save(groupHandle, "vertices", graph.numberOfVertices());
        int ID = GraphTraitsHDF5<Graph>::ID;
        save(groupHandle, "graph-type-id", ID);
    } catch (std::exception& e) {
        sError = e.what();
    }
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "CompleteGraph: Saving to dataset '"+graphName+"' failed: "+sError
        );
    }
}

template <class VISITOR>
void
load(
    const hid_t parentHandle,
    const std::string& graphName,
    CompleteGraph<VISITOR>& graph
) {
    typedef CompleteGraph<VISITOR> CompleteGraph;
    hdf5::HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName);
    std::string sError;
    try{
        {
            int ID;
            load(groupHandle, "graph-type-id", ID);
            if(ID!=GraphTraitsHDF5<CompleteGraph>::ID) {
                sError = "Stored graph type is not a CompleteGraph.";
                goto cleanup;
            }
        }
        std::size_t numberOfVertices;
        load(groupHandle, "vertices", numberOfVertices);
        graph.assign(numberOfVertices);
        }catch(std::exception& e) {
            sError = e.what();
    }
cleanup:
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "CompleteGraph: Loading dataset '"+graphName+"' failed: "+sError
        );
    }
}


} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_COMPLETE_GRAPH_HDF5_HXX