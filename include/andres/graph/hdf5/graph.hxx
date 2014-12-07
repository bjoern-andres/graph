#pragma once
#ifndef ANDRES_GRAPH_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_GRAPH_HDF5_HXX

#include <array>

#include "hdf5.hxx"

#include "../graph.hxx"
#include "mutable-graph.hxx"


namespace andres{
namespace graph{
namespace hdf5{

template<>
template<class VISITOR>
struct graph_traits<Graph<VISITOR> > {
    static const int ID = 10000;
};

template <class VISITOR>
void save(const hid_t, const std::string&, const Graph<VISITOR>&);

template <class VISITOR>
void load(const hid_t, const std::string&, Graph<VISITOR>&);

// Implementation

template <class VISITOR>
void
save(
    const hid_t fileHandle,
    const std::string& datasetName,
    const Graph<VISITOR>& graph
) {
    saveMutableGraph(fileHandle, datasetName, graph);
}

template <class VISITOR>
void
load(
    const hid_t fileHandle,
    const std::string& datasetName,
    Graph<VISITOR>& graph
) {
    loadMutableGraph(fileHandle, datasetName, graph);
}

} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_GRAPH_HDF5_HXX