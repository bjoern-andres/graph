#pragma once
#ifndef ANDRES_GRAPH_HDF5_MUTABLE_GRAPH
#define ANDRES_GRAPH_HDF5_MUTABLE_GRAPH

#include "hdf5.hxx"

namespace andres{
namespace graph{
namespace hdf5{

struct DatasetNamesMutable: public DatasetNames {
    const char vertices[9] = "vertices";
    const char edgeIndices[13] = "edge-indices";
    const char multipleEdgesEnabled[23] = "multiple-edges-enabled";
};
static const DatasetNamesMutable datasetNamesMutable;

template <class GRAPH>
void saveMutableGraph(const hid_t, const std::string&, const GRAPH&);

template <class GRAPH>
void loadMutableGraph(const hid_t, const std::string&, GRAPH&);

// Implementation

template <class GRAPH>
void
saveMutableGraph(
    const hid_t fileHandle,
    const std::string& datasetName,
    const GRAPH& graph
) {
    typedef GRAPH Graph;
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(fileHandle, datasetName,true);
    
    try {
        {
            int ID = graph_traits<Graph>::ID;
            save(groupHandle, datasetNamesMutable.graphTypeId, ID);
        }
        {
            const unsigned char multipleEdgesEnabled = graph.multipleEdgesEnabled();
            save(groupHandle, datasetNamesMutable.multipleEdgesEnabled, multipleEdgesEnabled);
        }
        
        const std::size_t numberOfEdges = graph.numberOfEdges();
        save(groupHandle, datasetNamesMutable.vertices, graph.numberOfVertices());
        
        std::vector<std::size_t> vecIJ;
        vecIJ.resize(2*numberOfEdges);
        std::size_t *ptrI = &vecIJ[0];
        std::size_t *ptrJ = &vecIJ[numberOfEdges];
        for(std::size_t e=0;e<numberOfEdges;++e) {
            *(ptrI++) = graph.vertexOfEdge(e,0);
            *(ptrJ++) = graph.vertexOfEdge(e,1);
        }
        {
            std::vector<std::size_t> shape = { numberOfEdges,2 };
            save(groupHandle, datasetNamesMutable.edgeIndices, shape, &vecIJ[0]);
        }
    } catch (std::exception& e) {
        closeGroup(groupHandle);
        throw std::runtime_error(
            "(Di)Graph: Save to dataset '"+datasetName+"'failed: " + std::string(e.what())
        );
    }
    closeGroup(groupHandle);
}

template <class GRAPH>
void
loadMutableGraph(
    const hid_t fileHandle,
    const std::string& datasetName,
    GRAPH& graph
) {
    typedef GRAPH Graph;
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(fileHandle, datasetName);
    
    std::string sError;
    
    try {
        {
            int ID;
            load(groupHandle, datasetNamesMutable.graphTypeId, ID);
            if(ID!=graph_traits<Graph>::ID) {
                sError = "Stored graph type mismatch.";
                goto cleanup;
            }
        }
        // Retrieve number of vertices
        std::size_t numberOfVertices;
        load(groupHandle, datasetNamesMutable.vertices, numberOfVertices);
        
        // Retrieve edge indices
        std::size_t numberOfEdges = 0;
        std::vector<std::size_t> bufferIJ;
        {
            std::vector<std::size_t> shape;
            load(groupHandle, datasetNamesMutable.edgeIndices, shape, bufferIJ);
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
            load(groupHandle, datasetNamesMutable.multipleEdgesEnabled, multipleEdgesEnabled);        
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
            "(Di)Graph: Load from dataset '"+datasetName+"' failed: " + sError
        );
    }
}

} //namespace hdf
} //namespace graph
} //namespace andres

#endif // ifndef ANDRES_GRAPH_HDF5_MUTABLE_GRAPH
