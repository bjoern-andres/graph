#pragma once
#ifndef ANDRES_GRAPH_CUT_VERTICES_HXX
#define ANDRES_GRAPH_CUT_VERTICES_HXX

#include <stack>
#include <vector>
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/subgraph.hxx"
#include "andres/graph/detail/cut-vertices.hxx"


namespace andres {
namespace graph {

template<class GRAPH>
struct CutVerticesBuffers {
    typedef GRAPH GraphType;

    CutVerticesBuffers(const GraphType&);

    std::vector<std::size_t> depth_;
    std::vector<char> is_cut_vertex_;
    std::vector<std::size_t> min_successor_depth_;
    std::vector<typename GraphType::VertexIterator> next_out_arc_;
    std::vector<int> parent_;
    std::vector<char> visited_;
};

template<typename GRAPH>
class CutVertices
{
public:
    typedef GRAPH GraphType;
    typedef CutVerticesBuffers<GraphType> CutVerticesBuffersType;

    
    CutVertices(const GraphType&);
    CutVertices(const GraphType&, CutVerticesBuffersType&);
    ~CutVertices();
    

    bool isCutVertex(std::size_t) const;
    
    void run();

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&);

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&, std::size_t);

private:
    CutVerticesBuffersType* buffer_;
    bool delete_buffer_;
    const GraphType& graph_;
};

template<typename GraphVisitor>
class CutVertices<CompleteGraph<GraphVisitor>>
{
public:
    typedef CompleteGraph<GraphVisitor> GraphType;
    typedef CutVerticesBuffers<GraphType> CutVerticesBuffersType;

    
    CutVertices(const GraphType&);
    CutVertices(const GraphType&, CutVerticesBuffersType&);
    ~CutVertices();
    

    bool isCutVertex(std::size_t) const;
    
    void run();

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&);

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&, std::size_t);

private:
    CutVerticesBuffersType* buffer_;
    bool delete_buffer_;
    const GraphType& graph_;
};




template<typename GRAPH>
CutVerticesBuffers<GRAPH>::CutVerticesBuffers(const GraphType& graph) :
    depth_(graph.numberOfVertices()),
    is_cut_vertex_(graph.numberOfVertices()),
    min_successor_depth_(graph.numberOfVertices()),
    next_out_arc_(graph.numberOfVertices()),
    parent_(graph.numberOfVertices()),
    visited_(graph.numberOfVertices())
{}





template<typename GRAPH>
inline
CutVertices<GRAPH>::CutVertices(const GraphType& graph) :
    buffer_(new CutVerticesBuffersType(graph)),
    delete_buffer_(true),
    graph_(graph)
{}

template<typename GRAPH>
inline
CutVertices<GRAPH>::CutVertices(const GraphType& graph, CutVerticesBuffersType& buffer) :  
    buffer_(&buffer),
    delete_buffer_(false),
    graph_(graph)
{}

template<typename GRAPH>
inline
CutVertices<GRAPH>::~CutVertices()
{
    if(delete_buffer_)
        delete buffer_;
}

template<typename GRAPH>
inline
bool CutVertices<GRAPH>::isCutVertex(std::size_t vertex_index) const
{
    return static_cast<bool>(buffer_->is_cut_vertex_[vertex_index]);
}

template<typename GRAPH>
inline
void CutVertices<GRAPH>::run()
{
    run(DefaultSubgraphMask<>());
}

template<typename GRAPH>
template<typename SUBGRAPH_MASK>
inline
void CutVertices<GRAPH>::run(const SUBGRAPH_MASK& subgraph_mask)
{
    std::fill(buffer_->parent_.begin(), buffer_->parent_.end(), -2);

    for (std::size_t i = 0; i < graph_.numberOfVertices(); ++i)
        if (buffer_->parent_[i] == -2 && subgraph_mask.vertex(i))
            run(subgraph_mask, i);
}

template<typename GRAPH>
template<typename SUBGRAPH_MASK>
inline
void CutVertices<GRAPH>::run(const SUBGRAPH_MASK& subgraph_mask, std::size_t starting_vertex)
{
    detail::find_cut_vertices(graph_, subgraph_mask, starting_vertex, *buffer_);
}







template<typename GraphVisitor>
inline
CutVertices<CompleteGraph<GraphVisitor>>::CutVertices(const GraphType& graph) :
    buffer_(new CutVerticesBuffersType(graph)),
    delete_buffer_(true),
    graph_(graph)
{}

template<typename GraphVisitor>
inline
CutVertices<CompleteGraph<GraphVisitor>>::CutVertices(const GraphType& graph, CutVerticesBuffersType& buffer) :  
    buffer_(&buffer),
    delete_buffer_(false),
    graph_(graph)
{}

template<typename GraphVisitor>
inline
CutVertices<CompleteGraph<GraphVisitor>>::~CutVertices()
{
    if(delete_buffer_)
        delete buffer_;
}

template<typename GraphVisitor>
inline
bool CutVertices<CompleteGraph<GraphVisitor>>::isCutVertex(std::size_t vertex_index) const
{
    return static_cast<bool>(buffer_->is_cut_vertex_[vertex_index]);
}

template<typename GraphVisitor>
inline
void CutVertices<CompleteGraph<GraphVisitor>>::run()
{
    std::fill(buffer_->is_cut_vertex_.begin(), buffer_->is_cut_vertex_.end(), 0);
}

template<typename GraphVisitor>
template<typename SUBGRAPH_MASK>
inline
void CutVertices<CompleteGraph<GraphVisitor>>::run(const SUBGRAPH_MASK& subgraph_mask)
{
    std::fill(buffer_->parent_.begin(), buffer_->parent_.end(), -2);

    for (std::size_t i = 0; i < graph_.numberOfVertices(); ++i)
        if (buffer_->parent_[i] == -2 && subgraph_mask.vertex(i))
            run(subgraph_mask, i);
}

template<typename GraphVisitor>
template<typename SUBGRAPH_MASK>
inline
void CutVertices<CompleteGraph<GraphVisitor>>::run(const SUBGRAPH_MASK& subgraph_mask, std::size_t starting_vertex)
{
    detail::find_cut_vertices(graph_, subgraph_mask, starting_vertex, *buffer_);
}


}
}

#endif