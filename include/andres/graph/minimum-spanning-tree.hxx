#pragma once
#ifndef ANDRES_GRAPH_MINIMUM_SPANNING_TREE_HXX
#define ANDRES_GRAPH_MINIMUM_SPANNING_TREE_HXX

#include <queue>
#include <stdexcept>
#include <vector>
#include "subgraph.hxx"
#include "detail/do_nothing_functor.hxx"


namespace andres {
namespace graph {

// Prim's algorithm to find a minimum spanning tree (forest) in an undirected graph
// 
// Robert C. Prim. (1957). Shortest connection networks and some generalizations
// Bell System Technical Journal, 36, pp. 1389–1401
// 
// Algorithm for sparse graphs has O(|E|log|V|) running time
// 
// A dynamic programming version for dense graphs has O(|V|^2) running time
//
template<typename GRAPH, typename ECA, typename PRED, typename FUNC = detail::do_nothing<typename ECA::value_type>>
inline 
typename ECA::value_type findMSTSparseGraph(const GRAPH& graph, const ECA& edge_weights, PRED& predecessor, const FUNC& f = FUNC());

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = detail::do_nothing<typename ECA::value_type>>
inline
typename ECA::value_type findMSTSparseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, PRED& predecessor, const FUNC& f = FUNC());

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = detail::do_nothing<typename ECA::value_type>>
inline
typename ECA::value_type findMSTSparseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, std::size_t starting_vertex, PRED& predecessor, const FUNC& f = FUNC());

template<typename GRAPH, typename ECA, typename PRED, typename FUNC = detail::do_nothing<typename ECA::value_type>>
inline 
typename ECA::value_type findMSTDenseGraph(const GRAPH& graph, const ECA& edge_weights, PRED& predecessor, const FUNC& f = FUNC());

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = detail::do_nothing<typename ECA::value_type>>
inline
typename ECA::value_type findMSTDenseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, PRED& predecessor, const FUNC& f = FUNC());

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = detail::do_nothing<typename ECA::value_type>>
inline
typename ECA::value_type findMSTDenseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, std::size_t starting_vertex, PRED& predecessor, const FUNC& f = FUNC());






template<typename GRAPH, typename ECA, typename PRED, typename FUNC>
inline 
typename ECA::value_type findMSTSparseGraph(const GRAPH& graph, const ECA& edge_weights, PRED& predecessor, const FUNC& f)
{
    return findMSTSparseGraph(graph, edge_weights, DefaultSubgraphMask<>(), predecessor, f);
}

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline
typename ECA::value_type findMSTSparseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, PRED& predecessor, const FUNC& f)
{
    typedef typename ECA::value_type value_type;

    std::fill(predecessor.begin(), predecessor.end(), graph.numberOfEdges());

    value_type mst_value = value_type();
    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (predecessor[i] == graph.numberOfEdges() && subgraph_mask.vertex(i))
            mst_value += findMSTSparseGraph(graph, edge_weights, subgraph_mask, i, predecessor, f);

    return mst_value;
}

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline
typename ECA::value_type findMSTSparseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, std::size_t starting_vertex, PRED& predecessor, const FUNC& f)
{
    typedef typename ECA::value_type value_type;

    struct element
    {
        element(value_type weight, std::size_t vertex_index) :
            weight_(weight), vertex_index_(vertex_index)
        {}

        bool operator<(const element& other) const
        {
            // as std::priority_queue is a max heap, invert comparison's logic
            return weight_ > other.weight_;
        }

        value_type weight_;
        std::size_t vertex_index_;
    };

    std::vector<value_type> min_edge(graph.numberOfVertices(), std::numeric_limits<value_type>::max());
    std::vector<char> visited(graph.numberOfVertices());

    std::priority_queue<element> Q;
    Q.push(element(value_type(), starting_vertex));

    predecessor[starting_vertex] = graph.numberOfEdges();

    value_type mst_value = value_type();

    while (!Q.empty())
    {
        auto v = Q.top();
        Q.pop();

        if (visited[v.vertex_index_])
            continue;

        visited[v.vertex_index_] = 1;
        mst_value += v.weight_;

        auto e = graph.edgesFromVertexBegin(v.vertex_index_);
        for (auto w = graph.verticesFromVertexBegin(v.vertex_index_); w != graph.verticesFromVertexEnd(v.vertex_index_); ++w, ++e)
            if (subgraph_mask.vertex(*w) &&
                subgraph_mask.edge(*e) &&
                !visited[*w] &&
                *w != v.vertex_index_ &&
                f(edge_weights[*e]) < min_edge[*w]
                )
            {
                min_edge[*w] = f(edge_weights[*e]);
                predecessor[*w] = *e;
                Q.push(element(f(edge_weights[*e]), *w));
            }
    }

    return mst_value;
}





template<typename GRAPH, typename ECA, typename PRED, typename FUNC>
inline 
typename ECA::value_type findMSTDenseGraph(const GRAPH& graph, const ECA& edge_weights, PRED& predecessor, const FUNC& f)
{
    return findMSTDenseGraph(graph, edge_weights, DefaultSubgraphMask<>(), predecessor, f);
}

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline
typename ECA::value_type findMSTDenseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, PRED& predecessor, const FUNC& f)
{
    typedef typename ECA::value_type value_type;

    std::fill(predecessor.begin(), predecessor.end(), graph.numberOfEdges());

    value_type mst_value = value_type();
    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (predecessor[i] == graph.numberOfEdges() && subgraph_mask.vertex(i))
            mst_value += findMSTDenseGraph(graph, edge_weights, subgraph_mask, i, predecessor, f);

    return mst_value;
}

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline
typename ECA::value_type findMSTDenseGraph(const GRAPH& graph, const ECA& edge_weights, const SUBGRAPH& subgraph_mask, std::size_t starting_vertex, PRED& predecessor, const FUNC& f)
{
    typedef typename ECA::value_type value_type;

    std::vector<value_type> min_edge(graph.numberOfVertices(), std::numeric_limits<value_type>::max());
    std::vector<char> visited(graph.numberOfVertices());

    min_edge[starting_vertex] = value_type();
    predecessor[starting_vertex] = graph.numberOfEdges();

    value_type mst_value = value_type();

    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
    {
        int v = -1;

        for (std::size_t j = 0; j < graph.numberOfVertices(); ++j)
            if (subgraph_mask.vertex(j) &&
                !visited[j] &&
                (v == -1 || min_edge[j] < min_edge[v])
                )
                v = j;

        if (v == -1)
            return mst_value;

        mst_value += min_edge[v];
        visited[v] = 1;

        auto e = graph.edgesFromVertexBegin(v);
        for (auto w = graph.verticesFromVertexBegin(v); w != graph.verticesFromVertexEnd(v); ++w, ++e)
            if (subgraph_mask.vertex(*w) &&
                subgraph_mask.edge(*e) &&
                !visited[*w] &&
                *w != v &&
                f(edge_weights[*e]) < min_edge[*w]
                )
            {
                min_edge[*w] = f(edge_weights[*e]);
                predecessor[*w] = *e;
            }
    }

    return mst_value;
}

}
}
#endif
