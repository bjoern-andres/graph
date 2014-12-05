#pragma once
#ifndef ANDRES_GRAPH_CUT_VERTICES_HXX
#define ANDRES_GRAPH_CUT_VERTICES_HXX

#include <stack>
#include <vector>
#include "andres/graph/subgraph.hxx"


namespace andres {
namespace graph {

// Hopcroft-Tarjan Algorithm to find cut vertices in an undirected graph
// 
// Hopcroft J. and Tarjan R. (1973). Efficient algorithms for graph manipulation
// Communications of the ACM 16 (6): 372–378
// 
template<class GRAPH>
struct CutVerticesBuffers
{
    typedef GRAPH GraphType;

    CutVerticesBuffers(const GraphType&);

    std::vector<std::size_t> depth_;
    std::vector<std::size_t> min_successor_depth_;
    std::vector<typename GraphType::VertexIterator> next_out_arc_;
    std::vector<int> parent_;
    std::vector<char> visited_;
};

template<typename GRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    std::vector<char>& isCutVertex
);

template<typename GRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    std::vector<char>& isCutVertex,
    CutVerticesBuffers<GRAPH>& buffer
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isCutVertex
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isCutVertex,
    CutVerticesBuffers<GRAPH>& buffer
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isCutVertex
);

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isCutVertex,
    CutVerticesBuffers<GRAPH>& buffer
);





template<typename GRAPH>
CutVerticesBuffers<GRAPH>::CutVerticesBuffers(const GraphType& graph) :
    depth_(graph.numberOfVertices()),
    min_successor_depth_(graph.numberOfVertices()),
    next_out_arc_(graph.numberOfVertices()),
    parent_(graph.numberOfVertices()),
    visited_(graph.numberOfVertices())
{}

template<typename GRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    std::vector<char>& isCutVertex
)
{
    auto buffer = CutVerticesBuffers<GRAPH>(graph);
    findCutVertices(graph, DefaultSubgraphMask<>(), isCutVertex, buffer);
}

template<typename GRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    std::vector<char>& isCutVertex,
    CutVerticesBuffers<GRAPH>& buffer
)
{
    findCutVertices(graph, DefaultSubgraphMask<>(), isCutVertex, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isCutVertex
)
{
    auto buffer = CutVerticesBuffers<GRAPH>(graph);
    findCutVertices(graph, subgraph_mask, isCutVertex, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::vector<char>& isCutVertex,
    CutVerticesBuffers<GRAPH>& buffer
)
{
    std::fill(buffer.parent_.begin(), buffer.parent_.end(), -2);

    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (buffer.parent_[i] == -2 && subgraph_mask.vertex(i))
            findCutVertices(graph, subgraph_mask, i, isCutVertex, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isCutVertex
)
{
    auto buffer = CutVerticesBuffers<GRAPH>(graph);
    findCutVertices(graph, subgraph_mask, starting_vertex, isCutVertex, buffer);
}

template<typename GRAPH, typename SUBGRAPH>
inline void
findCutVertices(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    std::vector<char>& isCutVertex,
    CutVerticesBuffers<GRAPH>& buffer
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
            isCutVertex[v] = 0;
            buffer.next_out_arc_[v] = graph.verticesFromVertexBegin(v);
            buffer.min_successor_depth_[v] = buffer.depth_[v];
        }
        else
        {
            auto to = *buffer.next_out_arc_[v];

            if (buffer.min_successor_depth_[to] >= buffer.depth_[v] && v != starting_vertex)
                isCutVertex[v] = 1;

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
                break;
            }
        }
    }

    std::size_t root_child_count = 0;
    for (auto to = graph.verticesFromVertexBegin(starting_vertex); to != graph.verticesFromVertexEnd(starting_vertex); ++to)
        if(buffer.parent_[*to] == starting_vertex)
            ++root_child_count;
        
    if(root_child_count >= 2)
        isCutVertex[starting_vertex] = 1;
}

}
}
#endif
