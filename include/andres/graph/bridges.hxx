#pragma once
#ifndef ANDRES_GRAPH_BRIDGES_HXX
#define ANDRES_GRAPH_BRIDGES_HXX

#include <stack>
#include <vector>
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/subgraph.hxx"
#include "andres/graph/detail/bridges.hxx"


namespace andres {
namespace graph {

template<class GRAPH>
struct BridgesBuffers {
    typedef GRAPH GraphType;

    BridgesBuffers(const GraphType&);

    std::vector<std::size_t> depth_;
    std::vector<char> is_bridge_;
    std::vector<std::size_t> min_successor_depth_;
    std::vector<typename GraphType::VertexIterator> next_out_arc_;
    std::vector<int> parent_;
    std::vector<char> visited_;
};

template<typename GRAPH>
class Bridges
{
public:
    typedef GRAPH GraphType;
    typedef BridgesBuffers<GraphType> BridgesBuffersType;

    Bridges(const GraphType&);
    Bridges(const GraphType&, BridgesBuffersType&);
    ~Bridges();

    bool isBridge(std::size_t) const;

    void run();

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&);

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&, std::size_t);

private:
    BridgesBuffersType* buffer_;
    bool delete_buffer_;
    const GraphType& graph_;
};

template<typename GraphVisitor>
class Bridges<CompleteGraph<GraphVisitor>>
{
public:
    typedef CompleteGraph<GraphVisitor> GraphType;
    typedef BridgesBuffers<GraphType> BridgesBuffersType;

    Bridges(const GraphType&);
    Bridges(const GraphType&, BridgesBuffersType&);
    ~Bridges();

    bool isBridge(std::size_t) const;

    void run();

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&);

    template<typename SUBGRAPH_MASK>
    void run(const SUBGRAPH_MASK&, std::size_t);

private:
    BridgesBuffersType* buffer_;
    bool delete_buffer_;
    const GraphType& graph_;
};





template<typename GRAPH>
BridgesBuffers<GRAPH>::BridgesBuffers(const GraphType& graph) :
    depth_(graph.numberOfVertices()),
    is_bridge_(graph.numberOfEdges()),
    min_successor_depth_(graph.numberOfVertices()),
    next_out_arc_(graph.numberOfVertices()),
    parent_(graph.numberOfVertices()),
    visited_(graph.numberOfVertices())
{}





template<typename GRAPH>
inline
Bridges<GRAPH>::Bridges(const GraphType& graph) :
    buffer_(new BridgesBuffersType(graph)),
    delete_buffer_(true),
    graph_(graph)
{}

template<typename GRAPH>
inline
Bridges<GRAPH>::Bridges(const GraphType& graph, BridgesBuffersType& buffer) :  
    buffer_(&buffer),
    delete_buffer_(false),
    graph_(graph)
{}

template<typename GRAPH>
inline
Bridges<GRAPH>::~Bridges()
{
    if(delete_buffer_)
        delete buffer_;
}

template<typename GRAPH>
inline
bool Bridges<GRAPH>::isBridge(std::size_t edge_index) const
{
    return static_cast<bool>(buffer_->is_bridge_[edge_index]);
}

template<typename GRAPH>
inline
void Bridges<GRAPH>::run()
{
    run(DefaultSubgraphMask<>());
}

template<typename GRAPH>
template<typename SUBGRAPH_MASK>
inline
void Bridges<GRAPH>::run(const SUBGRAPH_MASK& subgraph_mask)
{
    std::fill(buffer_->parent_.begin(), buffer_->parent_.end(), -2);

    for (std::size_t i = 0; i < graph_.numberOfVertices(); ++i)
        if (buffer_->parent_[i] == -2 && subgraph_mask.vertex(i))
            run(subgraph_mask, i);
}

template<typename GRAPH>
template<typename SUBGRAPH_MASK>
inline
void Bridges<GRAPH>::run(const SUBGRAPH_MASK& subgraph_mask, std::size_t starting_vertex)
{
    detail::find_bridges(graph_, subgraph_mask, starting_vertex, *buffer_);
}







template<typename GraphVisitor>
inline
Bridges<CompleteGraph<GraphVisitor>>::Bridges(const GraphType& graph) :
    buffer_(new BridgesBuffersType(graph)),
    delete_buffer_(true),
    graph_(graph)
{}

template<typename GraphVisitor>
inline
Bridges<CompleteGraph<GraphVisitor>>::Bridges(const GraphType& graph, BridgesBuffersType& buffer) :  
    buffer_(&buffer),
    delete_buffer_(false),
    graph_(graph)
{}

template<typename GraphVisitor>
inline
Bridges<CompleteGraph<GraphVisitor>>::~Bridges()
{
    if(delete_buffer_)
        delete buffer_;
}

template<typename GraphVisitor>
inline
bool Bridges<CompleteGraph<GraphVisitor>>::isBridge(std::size_t edge_index) const
{
    return static_cast<bool>(buffer_->is_bridge_[edge_index]);
}
template<typename GraphVisitor>
inline
void Bridges<CompleteGraph<GraphVisitor>>::run()
{
    std::fill(buffer_->is_bridge_.begin(), buffer_->is_bridge_.end(), 0);
}

template<typename GraphVisitor>
template<typename SUBGRAPH_MASK>
inline
void Bridges<CompleteGraph<GraphVisitor>>::run(const SUBGRAPH_MASK& subgraph_mask)
{
    std::fill(buffer_->parent_.begin(), buffer_->parent_.end(), -2);

    for (std::size_t i = 0; i < graph_.numberOfVertices(); ++i)
        if (buffer_->parent_[i] == -2 && subgraph_mask.vertex(i))
            run(subgraph_mask, i);
}

template<typename GraphVisitor>
template<typename SUBGRAPH_MASK>
inline
void Bridges<CompleteGraph<GraphVisitor>>::run(const SUBGRAPH_MASK& subgraph_mask, std::size_t starting_vertex)
{
    detail::find_bridges(graph_, subgraph_mask, starting_vertex, *buffer_);
}


}
}

#endif