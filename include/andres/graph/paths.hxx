#pragma once
#ifndef ANDRES_GRAPH_PATHS_HXX
#define ANDRES_GRAPH_PATHS_HXX

#include <cstddef>
#include <utility> // std::pair

#include "andres/graph/graph.hxx" // DefaultSubgraphMask

namespace andres {
namespace graph {

/// Search a path for a chord.
///
/// \param graph Graph.
/// \param begin Iterator to the beginning of the sequence of nodes on the path.
/// \param end Iterator to the end of the sequence of nodes on the path.
/// \param ignoreEdgeBetweenFirstAndLast Flag.
///
template<class GRAPH, class ITERATOR>
inline std::pair<bool, std::size_t>
findChord(
    const GRAPH& graph,
    ITERATOR begin,
    ITERATOR end,
    const bool ignoreEdgeBetweenFirstAndLast = false
) {
    return findChord(graph, DefaultSubgraphMask<>(), begin, end, 
        ignoreEdgeBetweenFirstAndLast);
}

/// Search a path for a chord.
///
/// \param graph Graph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param begin Iterator to the beginning of the sequence of nodes on the path.
/// \param end Iterator to the end of the sequence of nodes on the path.
/// \param ignoreEdgeBetweenFirstAndLast Flag.
///
template<class GRAPH, class SUBGRAPH_MASK, class ITERATOR>
inline std::pair<bool, std::size_t>
findChord(
    const GRAPH& graph,
    const SUBGRAPH_MASK& mask,
    ITERATOR begin,
    ITERATOR end,
    const bool ignoreEdgeBetweenFirstAndLast = false
) {
    for(ITERATOR it = begin; it != end - 1; ++it) 
    for(ITERATOR it2 = it + 2; it2 != end; ++it2) {
        if(ignoreEdgeBetweenFirstAndLast && it == begin && it2 == end - 1) {
            continue;
        }
        std::pair<bool, std::size_t> p = graph.findEdge(*it, *it2);
        if(p.first && mask.edge(p.second)) {
            return p;
        }
    }
    return std::pair<bool, std::size_t>(false, 0);
}

/// Determine whether a path has a chord in time linear in the path length.
///
/// \param graph Graph.
/// \param begin Iterator to the beginning of the sequence of nodes on the path.
/// \param end Iterator to the end of the sequence of nodes on the path.
/// \param seen Vector of labels (should contain initially zeros)
/// \param ignoreEdgeBetweenFirstAndLast Flag.
///
template<class GRAPH, class ITERATOR>
inline bool 
hasChord(
    const GRAPH& graph,
    ITERATOR begin,
    ITERATOR end,
    std::vector<char>& seen,
    const bool ignoreEdgeBetweenFirstAndLast = false)
{
    for (ITERATOR path_it = begin; path_it != end; path_it++)
    {
        if (seen[*path_it])
            return true;

        for (auto adj_it = graph.adjacenciesFromVertexBegin(*path_it); adj_it != graph.adjacenciesFromVertexEnd(*path_it); adj_it++)
        {
            // exclude edge between first and last vertex in the path
            if (ignoreEdgeBetweenFirstAndLast && path_it == begin && adj_it->vertex() == *(end-1))
                continue;

            if (adj_it->vertex() == *(path_it+1) )
                continue;
            
            seen[adj_it->vertex()] = 1;
        }
    }
    return false;
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_PATHS_HXX
