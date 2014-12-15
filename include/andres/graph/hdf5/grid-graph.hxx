#pragma once
#ifndef ANDRES_GRAPH_GRID_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_GRID_GRAPH_HDF5_HXX

#include <stdexcept>
#include <string>

#include "hdf5.hxx"
#include "../grid-graph.hxx"


namespace andres{
namespace graph{
namespace hdf5{

template<>
template<unsigned char D, class VISITOR>
struct GraphTraitsHDF5<GridGraph<D, VISITOR> > {
    static const int ID = 10003;
};

template<unsigned char D, class VISITOR>
void save(const hid_t, const std::string&, const GridGraph<D, VISITOR>&);

template<unsigned char D, class VISITOR>
void load(const hid_t, const std::string&, GridGraph<D, VISITOR>&);

// Implementation

template<unsigned char D, class VISITOR>
void
save(
    const hid_t parentHandle,
    const std::string& graphName,
    const GridGraph<D, VISITOR>& graph
) {
    typedef GridGraph<D, VISITOR> Graph;
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName,true);
    std::string sError;
    try {
        const int ID = GraphTraitsHDF5<Graph>::ID;
        save(groupHandle, "graph-type-id", ID);
        // Save shape
        std::array<std::size_t, D> shape;
        for(std::size_t i=0;i<D;++i) {
            shape[i] = graph.shape(i);
        }
        save(groupHandle, "shape", {D}, &shape[0]);
    } catch(std::exception& e) {
        sError = e.what();
    }
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "GridGraph: Save to dataset '"+graphName+"' failed :"+sError
        );
    }
}

template<unsigned char D, class VISITOR>
void
load(
    const hid_t parentHandle,
    const std::string& graphName,
    GridGraph<D, VISITOR>& graph
) {
    typedef GridGraph<D, VISITOR> Graph;
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName);
    
    std::string sError;
    try{
        std::vector<std::size_t> nDims;
        std::vector<std::size_t> shape;
        {
            int ID;
            load(groupHandle, "graph-type-id", ID);
            if(ID != GraphTraitsHDF5<Graph>::ID) {
                sError = "Stored graph type is not a GridGraph.";
                goto cleanup;
            }
        }
        load(groupHandle, "shape", nDims, shape);
        if(nDims.size() != 1) {
            sError = "Shape dataset is not a vector.";
            goto cleanup;
        }
        if(nDims[0] != D) {
            sError = "Saved Graph has a mismatching dimension.";
            goto cleanup;
        }
        typename GridGraph<D, VISITOR>::VertexCoordinate vCoord;
        std::copy(shape.begin(), shape.end(), vCoord.begin());
        graph.assign(vCoord);
    } catch(std::exception& e) {
        sError = "Exception during loading: " + std::string(e.what());
        goto cleanup;
    }
cleanup:
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "GridGraph: Loading from dataset '" +graphName+"' failed: "+sError
        );
    }
}


} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_GRID_GRAPH_HDF5_HXX
