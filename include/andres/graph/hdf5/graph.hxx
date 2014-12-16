#pragma once
#ifndef ANDRES_GRAPH_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_GRAPH_HDF5_HXX

#include <array>

#include "hdf5.hxx"

#include "../graph.hxx"


namespace andres{
namespace graph{
namespace hdf5{

template<>
template<class VISITOR>
struct GraphTraitsHDF5<Graph<VISITOR> > {
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
    const hid_t parentHandle,
    const std::string& graphName,
    const Graph<VISITOR>& graph
) {
    typedef Graph<VISITOR> Graph;
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName,true);
    
    try {
        {
            int ID = GraphTraitsHDF5<Graph>::ID;
            save(groupHandle, "graph-type-id", ID);
        }
        {
            const unsigned char multipleEdgesEnabled = graph.multipleEdgesEnabled();
            save(groupHandle, "multiple-edges-enabled", multipleEdgesEnabled);
        }
        
        const std::size_t numberOfEdges = graph.numberOfEdges();
        save(groupHandle, "vertices", graph.numberOfVertices());
        
        std::vector<std::size_t> vecIJ;
        vecIJ.resize(2*numberOfEdges);
        std::size_t *ptrI = &vecIJ[0];
        std::size_t *ptrJ = &vecIJ[numberOfEdges];
        for(std::size_t e=0;e<numberOfEdges;++e) {
            *(ptrI++) = graph.vertexOfEdge(e,0);
            *(ptrJ++) = graph.vertexOfEdge(e,1);
        }
        save(groupHandle, "edge-indices", {numberOfEdges, 2}, &vecIJ[0]);
    } catch (std::exception& e) {
        closeGroup(groupHandle);
        throw std::runtime_error(
            "(Di)Graph: Save to dataset '"+graphName+"'failed: " + std::string(e.what())
        );
    }
    closeGroup(groupHandle);
}

template <class VISITOR>
void
load(
    const hid_t parentHandle,
    const std::string& graphName,
    Graph<VISITOR>& graph
) {
    typedef Graph<VISITOR> Graph;
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName);
    
    std::string sError;
    
    try {
        {
            int ID;
            load(groupHandle, "graph-type-id", ID);
            if(ID!=GraphTraitsHDF5<Graph>::ID) {
                sError = "Stored graph type mismatch.";
                goto cleanup;
            }
        }
        // Retrieve number of vertices
        std::size_t numberOfVertices;
        load(groupHandle, "vertices", numberOfVertices);
        
        // Retrieve edge indices
        std::size_t numberOfEdges = 0;
        std::vector<std::size_t> bufferIJ;
        {
            std::vector<std::size_t> shape;
            load(groupHandle, "edge-indices", shape, bufferIJ);
            if(shape.size()!=2) {
                sError = "Stored edges are not 2-dimensional.";
                goto cleanup;
            }
            if(shape[1]!=2) {
                sError = "Stored number of edge endpoints is not 2.";
                goto cleanup;
            }
            assert(shape[0]*2 == bufferIJ.size());
            numberOfEdges = shape[0];
        }
        std::size_t* ptrI = &bufferIJ[0];
        std::size_t* ptrJ = &bufferIJ[numberOfEdges];
        graph.assign(numberOfVertices);
        {
            unsigned char multipleEdgesEnabled;
            load(groupHandle, "multiple-edges-enabled", multipleEdgesEnabled);        
            graph.multipleEdgesEnabled() = multipleEdgesEnabled;
        }
        graph.reserveEdges(numberOfEdges);
        {
            for(std::size_t i=0;i<numberOfEdges;++i) {
                const std::size_t s = *(ptrI++);
                const std::size_t t = *(ptrJ++);
                if(s>=numberOfVertices || t>=numberOfVertices) {
                    std::stringstream s;
                    s << "Vertex index of edge (" << i << ") out of bounds.";
                    sError = s.str();
                    goto cleanup;
                }
                const std::size_t e = graph.insertEdge(s,t);
                assert(e == graph.numberOfEdges()-1);
            }
        }
    } catch(std::exception& e) {
        sError = e.what();
    }
cleanup:
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error(
            "(Di)Graph: Load from dataset '"+graphName+"' failed: " + sError
        );
    }
}

} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_GRAPH_HDF5_HXX