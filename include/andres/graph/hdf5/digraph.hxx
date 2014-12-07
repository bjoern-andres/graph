#pragma once
#ifndef ANDRES_GRAPH_DIGRAPH_HDF5_HXX
#define ANDRES_GRAPH_DIGRAPH_HDF5_HXX

#include "hdf5.hxx"
#include "../digraph.hxx"
#include "mutable-graph.hxx"

namespace andres{
namespace graph{
namespace hdf5{

template<>
template<class VISITOR>
struct graph_traits<Digraph<VISITOR> > {
    static const int ID = 10001;
};

template <class VISITOR>
void save(const hid_t, const std::string&, const Digraph<VISITOR>&);

template <class VISITOR>
void load(const hid_t, const std::string&, Digraph<VISITOR>&);

// Implementation

template <class VISITOR>
void
save(
    const hid_t fileHandle,
    const std::string& datasetName,
    const Digraph<VISITOR>& graph
) {
    saveMutableGraph(fileHandle, datasetName, graph);
}

template <class VISITOR>
void
load(
    const hid_t fileHandle,
    const std::string& datasetName,
    Digraph<VISITOR>& graph
) {
    loadMutableGraph(fileHandle, datasetName, graph);
}


} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_DIGRAPH_HDF5_HXX