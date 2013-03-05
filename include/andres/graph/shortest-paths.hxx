// Copyright (c) 2013 by Bjoern Andres.
// 
// This software was developed by Bjoern Andres.
// Enquiries shall be directed to bjoern@andres.sc.
//
// All advertising materials mentioning features or use of this software must
// display the following acknowledgement: ``This product includes andres::graph
// developed by Bjoern Andres. Please direct enquiries concerning andres::graph
// to bjoern@andres.sc''.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// - All advertising materials mentioning features or use of this software must 
//   display the following acknowledgement: ``This product includes 
//   andres::graph developed by Bjoern Andres. Please direct enquiries 
//   concerning andres::graph to bjoern@andres.sc''.
// - The name of the author must not be used to endorse or promote products 
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#pragma once
#ifndef ANDRES_GRAPH_SORTEST_PATHS_HXX
#define ANDRES_GRAPH_SORTEST_PATHS_HXX

#include <limits> // std::numeric_limits
#include <deque>
#include <queue>
#include <vector>

#include "andres/graph/graph.hxx" // DefaultSubgraphMask

namespace andres {
namespace graph {

// \cond SUPPRESS_DOXYGEN
namespace graph_detail {

template<class T>
inline void
spspHelper(
    const std::vector<ptrdiff_t>& parents,
    const T vPositive,
    const T vNegative,
    std::deque<T>& path
) {
    assert(vPositive >= 0);
    assert(vNegative >= 0);
    T t = vPositive;
    for(;;) {
        path.push_front(t);
        if(parents[t] - 1 == t) {
            break;
        }
        else {
            t = parents[t] - 1;
        }
    }
    t = vNegative;
    for(;;) {
        path.push_back(t);
        if(-parents[t] - 1 == t) {
            break;
        }
        else {
            t = -parents[t] - 1;
        }
    }
}

template<class T>
struct DijkstraQueueEntry {
    typedef T Value;

    DijkstraQueueEntry(const size_t vertex = 0, const Value distance = Value())
        :   vertex_(vertex), distance_(distance)
        {}
    bool operator<(const DijkstraQueueEntry<Value>& other) const
        { return distance_ > other.distance_; }
    bool operator==(const DijkstraQueueEntry<Value>& other) const
        { return vertex_ == other.vertex_ && distance_ == other.distance_; }
    bool operator!=(const DijkstraQueueEntry<Value>& other) const
        { return vertex_ != other.vertex_ || distance_ != other.distance_; }

    size_t vertex_;
    Value distance_;
};

} // namespace graph_detail
// \endcond

/// Search for a shortest path connecting a single pair of vertices in an unweighted graph using breadth-first-search.
///
/// This function works for both undirected and directed graphs. It carries out
/// breadth-first searches from the source vertex vs and the target vertex vt, 
/// alternating between the two search trees, until these trees meet and thus, 
/// a shortest path from vs to vt has been found.
///
/// \param graph A graph class such as andres::Graph or andres::Digraph.
/// \param vs The source vertex.
/// \param vt The target vertex.
/// \param path A double-ended queue to which the path is written.
/// \param parents An optional external buffer.
/// \return true if a (shortest) path was found, false otherwise.
///
template<class GRAPH>
inline bool 
spsp(
    const GRAPH& graph, 
    const size_t vs,
    const size_t vt,
    std::deque<size_t>& path,
    std::vector<ptrdiff_t>& parents = std::vector<ptrdiff_t>()
) {
    return spsp(graph, DefaultSubgraphMask<>(), vs, vt, path, parents);
}

/// Search for a shortest path connecting a single pair of vertices in an unweighted subgraph using breadth-first-search.
///
/// This function works for both undirected and directed graphs. It carries out
/// breadth-first searches from the source vertex vs and the target vertex vt, 
/// alternating between the two search trees, until these trees meet and thus, 
/// a shortest path from vs to vt has been found.
///
/// \param graph A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs The source vertex.
/// \param vt The target vertex.
/// \param path A double-ended queue to which the path is written.
/// \param parents An optional external buffer.
/// \return true if a (shortest) path was found, false otherwise.
///
template<class GRAPH, class SUBGRAPH_MASK>
bool 
spsp(
    const GRAPH& graph, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const size_t vt,
    std::deque<size_t>& path,
    std::vector<ptrdiff_t>& parents = std::vector<ptrdiff_t>()
) {
    path.clear();
    if(!mask.vertex(vs) || !mask.vertex(vt)) {
        return false;
    }
    if(vs == vt) {
        path.push_back(vs);
        return true;
    }
    parents.resize(graph.numberOfVertices());
    std::fill(parents.begin(), parents.end(), 0);
    parents[vs] = vs + 1;
    parents[vt] = -static_cast<ptrdiff_t>(vt) - 1;
    std::queue<size_t> queues[2];
    queues[0].push(vs);
    queues[1].push(vt);
    for(size_t q = 0; true; q = 1 - q) { // infinite loop, alternating queues
        const size_t numberOfNodesAtFront = queues[q].size();
        for(size_t n = 0; n < numberOfNodesAtFront; ++n) {
            const size_t v = queues[q].front();
            queues[q].pop();
            GRAPH::AdjacencyIterator it;
            GRAPH::AdjacencyIterator end;
            if(q == 0) {
                it = graph.adjacenciesFromVertexBegin(v);
                end = graph.adjacenciesFromVertexEnd(v);
            }
            else {
                it = graph.adjacenciesToVertexBegin(v);
                end = graph.adjacenciesToVertexEnd(v);
            }
            for(; it != end; ++it) {
                if(!mask.edge(it->edge()) || !mask.vertex(it->vertex())) {
                    continue;
                }
                if(parents[it->vertex()] < 0 && q == 0) {
                    graph_detail::spspHelper(parents, v, it->vertex(), path);
                    assert(path[0] == vs);
                    assert(path.back() == vt);
                    return true;
                }
                else if(parents[it->vertex()] > 0 && q == 1) {
                    graph_detail::spspHelper(parents, it->vertex(), v, path);
                    assert(path[0] == vs);
                    assert(path.back() == vt);
                    return true;
                }
                else if(parents[it->vertex()] == 0) {
                    if(q == 0) {
                        parents[it->vertex()] = v + 1;
                    }
                    else {
                        parents[it->vertex()] = -static_cast<ptrdiff_t>(v) - 1;
                    }
                    queues[q].push(it->vertex());
                }
            }
        }
        if(queues[0].empty() && queues[1].empty()) {
            return false;
        }
    }
}

/// Search for shortest paths from a given vertex to every other vertex in a graph with unit edge weights using Dijkstra's algorithm.
///
/// \param graph A graph class such as andres::Graph or andres::Digraph.
/// \param vs Source vertex.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class DISTANCES, class PARENTS>
inline void
sssp(
    const GRAPH& graph, 
    const size_t vs,
    DISTANCES& distances,
    PARENTS& parents = std::vector<size_t>(graph.numberOfVertices())
) {
    typedef typename std::iterator_traits<DISTANCES>::value_type Value;
    sssp(graph, DefaultSubgraphMask<>(), vs,
        UnitEdgeWeightIterator<Value>(), distances, parents
    );
}

/// Search for shortest paths from a given vertex to every other vertex in a subgraph with unit edge weights using Dijkstra's algorithm.
///
/// \param graph A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class SUBGRAPH_MASK, class DISTANCES, class PARENTS>
inline void 
sssp(
    const GRAPH& graph, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    DISTANCES& distances,
    PARENTS& parents = std::vector<size_t>(graph.numberOfVertices())
) {
    typedef typename std::iterator_traits<DISTANCES>::value_type Value;
    sssp(graph, mask, vs, UnitEdgeWeightIterator<Value>(), 
        distances, parents
    );
}

/// Search for shortest paths from a given vertex to every other vertex in a graph with non-negative edge weights using Dijkstra's algorithm.
///
/// \param graph A graph class such as andres::Graph or andres::Digraph.
/// \param vs Source vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class EDGE_WEIGHTS, class DISTANCES, class PARENTS>
inline void
sssp(
    const GRAPH& graph, 
    const size_t vs,
    const EDGE_WEIGHTS& edgeWeights,
    DISTANCES& distances,
    PARENTS& parents = std::vector<size_t>(graph.numberOfVertices())
) {
    sssp(graph, DefaultSubgraphMask<>(), vs, edgeWeights, distances, parents);
}

/// Search for shortest paths from a given vertex to every other vertex in a subgraph with non-negative edge weights using Dijkstra's algorithm.
///
/// \param graph A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class SUBGRAPH_MASK, class EDGE_WEIGHTS, class DISTANCES, class PARENTS>
void 
sssp(
    const GRAPH& graph, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const EDGE_WEIGHTS& edgeWeights,
    DISTANCES& distances,
    PARENTS& parents = std::vector<size_t>(graph.numberOfVertices())
) {
    typedef std::iterator_traits<DISTANCES>::value_type Value;
    typedef graph_detail::DijkstraQueueEntry<Value> Entry;

    assert(mask.vertex(vs));  
    const Value infinity = std::numeric_limits<Value>::has_infinity 
        ? std::numeric_limits<Value>::infinity() 
        : std::numeric_limits<Value>::max();
    std::priority_queue<Entry> queue;
    distances[vs] = 0;
    queue.push(vs);
    for(size_t v = 0; v < graph.numberOfVertices(); ++v) {
        if(mask.vertex(v) && v != vs) {
            distances[v] = infinity;
            queue.push(Entry(v, infinity));
        }
    }
    while(!queue.empty()) {
        const size_t v = queue.top().vertex_;
        queue.pop();
        if(distances[v] == infinity) {
            return;
        }
        for(GRAPH::AdjacencyIterator it = graph.adjacenciesFromVertexBegin(v);
        it != graph.adjacenciesFromVertexEnd(v); ++it) {
            if(mask.vertex(it->vertex()) && mask.edge(it->edge())) {
                const Value alternativeDistance = distances[v] + edgeWeights[it->edge()];
                if(alternativeDistance < distances[it->vertex()]) {
                    distances[it->vertex()] = alternativeDistance;
                    parents[it->vertex()] = v;
                    queue.push(Entry(it->vertex(), alternativeDistance));
                    // pushing v another time, not worring about existing entries
                    // v in the queue at deprecated positions.
                    // alternatively, one could use a heap from which elements 
                    // can be removed. this is beneficial for dense graphs in 
                    // which many weights are equal.
                }
            }
        }
    }
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_SORTEST_PATHS_HXX
