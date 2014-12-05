#pragma once
#ifndef ANDRES_GRAPH_BRIDGES_HXX
#define ANDRES_GRAPH_BRIDGES_HXX

#include <stack>
#include <vector>
#include "andres/graph/subgraph.hxx"


namespace andres {
namespace graph {

// Tarjan's Algorithm to find bridges in an undirected graph
//
// Tarjan, R. (1974). A note on finding the bridges of a graph
// Information Processing Letters 2 (6): 160–161
// 
template<class GRAPH>
struct BridgesBuffers
{
    typedef GRAPH GraphType;

    BridgesBuffers(const GraphType&);

    std::vector<std::size_t> depth_;
    std::vector<std::size_t> min_successor_depth_;
    std::vector<typename GraphType::VertexIterator> next_out_arc_;
    std::vector<int> parent_;
    std::vector<char> visited_;
};

template<typename GRAPH>
inline void
findBridges(
    const GRAPH& graph,
    std::vector<char>& isBridge
);

template<typename GRAPH>
inline void
findBridges(
    const GRAPH& graph,
    std::vector<char>& isBridge,
    BridgesBuffers<GRAPH>& buffer
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isBridge
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isBridge,
    BridgesBuffers<GRAPH>& buffer
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isBridge
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isBridge,
    BridgesBuffers<GRAPH>& buffer
);






template<typename GRAPH>
BridgesBuffers<GRAPH>::BridgesBuffers(const GraphType& graph) :
    depth_(graph.numberOfVertices()),
    min_successor_depth_(graph.numberOfVertices()),
    next_out_arc_(graph.numberOfVertices()),
    parent_(graph.numberOfVertices()),
    visited_(graph.numberOfVertices())
{}

template<typename GRAPH>
inline void
findBridges(
    const GRAPH& graph,
    std::vector<char>& isBridge
)
{
    auto buffer = BridgesBuffers<GRAPH>(graph);
    findBridges(graph, DefaultSubgraphMask<>(), isBridge, buffer);
}

template<typename GRAPH>
inline void
findBridges(
    const GRAPH& graph,
    std::vector<char>& isBridge,
    BridgesBuffers<GRAPH>& buffer
)
{
    findBridges(graph, DefaultSubgraphMask<>(), isBridge, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isBridge
)
{
    auto buffer = BridgesBuffers<GRAPH>(graph);
    findBridges(graph, subgraph_mask, isBridge, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isBridge,
    BridgesBuffers<GRAPH>& buffer
)
{
   std::fill(buffer.parent_.begin(), buffer.parent_.end(), -2);

    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (buffer.parent_[i] == -2 && subgraph_mask.vertex(i))
            findBridges(graph, subgraph_mask, i, isBridge, buffer); 
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isBridge
)
{
    auto buffer = BridgesBuffers<GRAPH>(graph);
    findBridges(graph, subgraph_mask, starting_vertex, isBridge, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isBridge,
    BridgesBuffers<GRAPH>& buffer
)
{
    std::fill(buffer.visited_.begin(), buffer.visited_.end(), 0);

    std::stack<std::size_t> S;

    S.push(starting_vertex);
    buffer.depth_[starting_vertex] = 0;
    buffer.parent_[starting_vertex] = -1;

    while (!S.empty())
    {
        auto v = S.top();
        S.pop();

        if (!buffer.visited_[v])
        {
            buffer.visited_[v] = 1;
            buffer.next_out_arc_[v] = graph.verticesFromVertexBegin(v);
            buffer.min_successor_depth_[v] = buffer.depth_[v];
        }
        else
        {
            auto to = *buffer.next_out_arc_[v];

            if (buffer.min_successor_depth_[to] > buffer.depth_[v])
            {
                typename GRAPH::EdgeIterator e = graph.edgesFromVertexBegin(v) + (buffer.next_out_arc_[v] - graph.verticesFromVertexBegin(v));
                isBridge[*e] = 1;
            }

            buffer.min_successor_depth_[v] = std::min(buffer.min_successor_depth_[v], buffer.min_successor_depth_[to]);
            ++buffer.next_out_arc_[v];
        }

        while (buffer.next_out_arc_[v] != graph.verticesFromVertexEnd(v))
        {
            typename GRAPH::EdgeIterator e = graph.edgesFromVertexBegin(v) + (buffer.next_out_arc_[v] - graph.verticesFromVertexBegin(v));

            if (
                !subgraph_mask.vertex(*buffer.next_out_arc_[v]) ||
                !subgraph_mask.edge(*e)
                )
            {
                ++buffer.next_out_arc_[v];
                continue;
            }

            auto to = *buffer.next_out_arc_[v];
            if (buffer.visited_[to])
            {
                if(buffer.parent_[v] != to)
                    buffer.min_successor_depth_[v] = std::min(buffer.min_successor_depth_[v], buffer.depth_[to]);

                ++buffer.next_out_arc_[v];
            }
            else
            {
                S.push(v);
                S.push(to);
                buffer.parent_[to] = v;
                buffer.depth_[to] = buffer.depth_[v] + 1;
                isBridge[*e] = 0;
                break;
            }
        }
    }
}

}
}

#endif