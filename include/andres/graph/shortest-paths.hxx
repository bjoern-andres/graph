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
#include <algorithm> // std::reverse

#include "andres/graph/graph.hxx" // DefaultSubgraphMask


namespace andres {
namespace graph {

template<class T>
inline void
spspHelper(
    const std::vector<std::ptrdiff_t>&,
    const T,
    const T,
    std::deque<T>&
);

template<class GRAPH>
inline bool
spsp(
    const GRAPH&,
    const size_t,
    const size_t,
    std::deque<size_t>&,
    std::vector<std::ptrdiff_t>&
);

template<class GRAPH>
inline bool
spsp(
    const GRAPH&,
    const size_t,
    const size_t,
    std::deque<size_t>&
);

template<class GRAPH, class SUBGRAPH_MASK>
bool
spsp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const size_t,
    std::deque<size_t>&
);

template<class GRAPH, class SUBGRAPH_MASK>
bool
spsp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const size_t,
    std::deque<size_t>&,
    std::vector<std::ptrdiff_t>&
);

template<
    class GRAPH,
    class EDGE_WEIGHT_ITERATOR,
    class T
>
inline void
spsp(
    const GRAPH&,
    const size_t,
    const size_t,
    EDGE_WEIGHT_ITERATOR,
    std::deque<size_t>&,
    T&
);

template<
    class GRAPH,
    class SUBGRAPH_MASK,
    class EDGE_WEIGHT_ITERATOR,
    class T
>
inline void
spsp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const size_t,
    EDGE_WEIGHT_ITERATOR,
    std::deque<size_t>&,
    T&
);

template<class GRAPH, class DISTANCE_ITERATOR>
inline void
sssp(
    const GRAPH&,
    const size_t,
    DISTANCE_ITERATOR
);

template<class GRAPH, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
sssp(
    const GRAPH&,
    const size_t,
    DISTANCE_ITERATOR,
    PARENT_ITERATOR
);

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR>
inline void
sssp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    DISTANCE_ITERATOR
);

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
sssp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    DISTANCE_ITERATOR,
    PARENT_ITERATOR
);

template<class GRAPH, class EDGE_WEIGHT_ITERATOR, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
sssp(
    const GRAPH&,
    const size_t,
    const EDGE_WEIGHT_ITERATOR,
    DISTANCE_ITERATOR,
    PARENT_ITERATOR
);

template<
    class GRAPH,
    class SUBGRAPH_MASK,
    class EDGE_WEIGHT_ITERATOR,
    class DISTANCE_ITERATOR,
    class PARENT_ITERATOR
>
inline void
sssp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const EDGE_WEIGHT_ITERATOR,
    DISTANCE_ITERATOR,
    PARENT_ITERATOR
);

template<
    class GRAPH,
    class SUBGRAPH_MASK,
    class EDGE_WEIGHT_ITERATOR,
    class DISTANCE_ITERATOR,
    class PARENT_ITERATOR,
    class VISITOR
>
void
sssp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const EDGE_WEIGHT_ITERATOR,
    DISTANCE_ITERATOR,
    PARENT_ITERATOR,
    VISITOR&
);

template<
    class GRAPH,
    class SUBGRAPH_MASK,
    class EDGE_WEIGHT_ITERATOR,
    class T,
    class DISTANCE_ITERATOR,
    class PARENT_ITERATOR
>
inline void
spsp(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const size_t,
    EDGE_WEIGHT_ITERATOR,
    std::deque<size_t>&,
    T&,
    DISTANCE_ITERATOR,
    PARENT_ITERATOR
);

// \cond SUPPRESS_DOXYGEN
namespace graph_detail {

template<class T>
inline void
spspHelper(
    const std::vector<std::ptrdiff_t>& parents,
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

// Single source shortest path visitor for Dijkstra's algorithm.
template<class DISTANCE_ITERATOR, class PARENT_ITERATOR>
class DijkstraSPSPVisitor {
public:
    typedef DISTANCE_ITERATOR Distances;
    typedef PARENT_ITERATOR Parents;
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;

    DijkstraSPSPVisitor(
        const size_t vs, 
        const size_t vt, 
        std::deque<size_t>& path
    )
    :   vs_(vs), vt_(vt), path_(path)
    {
        path_.clear(); // so the path will be empty if no path is found
    }
    
    bool operator()(
        Distances distances, 
        Parents parents, 
        size_t vertex
    ) {
        if(vertex == vt_) {
            const Value infinity = std::numeric_limits<Value>::has_infinity 
                ? std::numeric_limits<Value>::infinity() 
                : std::numeric_limits<Value>::max();
            if(distances[vertex] == infinity) {
                path_.clear();
                return false; // stop the algorithm
            }
            for(;;) {
                path_.push_front(vertex);
                if(vertex == vs_) {
                    return false; // stop the algorithm
                }
                else {
                    vertex = parents[vertex];
                }
            }
        }
        else {
            return true; // continue with the algorithm
        }
    }

    size_t vs_;
    size_t vt_;
    std::deque<size_t>& path_;
};

} // namespace graph_detail
// \endcond

/// Search for a shortest path from one to another vertex in an **unweighted** graph using breadth-first-search.
///
/// This function works for both undirected and directed graphs. It carries out
/// breadth-first searches from the source vertex vs and the target vertex vt, 
/// alternating between the two search trees, until these trees meet and thus, 
/// a shortest path from vs to vt has been found.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param vs The source vertex.
/// \param vt The target vertex.
/// \param path A double-ended queue to which the path is written.
/// \param parents An optional external buffer.
/// \return true if a (shortest) path was found, false otherwise.
///
template<class GRAPH>
inline bool 
spsp(
    const GRAPH& g, 
    const size_t vs,
    const size_t vt,
    std::deque<size_t>& path,
    std::vector<std::ptrdiff_t>& parents
) {
    return spsp(g, DefaultSubgraphMask<>(), vs, vt, path, parents);
}

template<class GRAPH>
inline bool
spsp(
    const GRAPH& g,
    const size_t vs,
    const size_t vt,
    std::deque<size_t>& path
) {
    std::vector<std::ptrdiff_t> parents = std::vector<std::ptrdiff_t>();
    return spsp(g, DefaultSubgraphMask<>(), vs, vt, path, parents);
}

/// Search for a shortest path from one to another vertex in an **unweighted subgraph** using breadth-first-search.
///
/// This function works for both undirected and directed graphs. It carries out
/// breadth-first searches from the source vertex vs and the target vertex vt, 
/// alternating between the two search trees, until these trees meet and thus, 
/// a shortest path from vs to vt has been found.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs The source vertex.
/// \param vt The target vertex.
/// \param path A double-ended queue to which the path is written.
/// \param parents An optional external buffer.
/// \return true if a (shortest) path was found, false otherwise.
///
///
///

template<class GRAPH, class SUBGRAPH_MASK>
bool
spsp(
    const GRAPH& g,
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const size_t vt,
    std::deque<size_t>& path
) {

    std::vector<std::ptrdiff_t> parents = std::vector<std::ptrdiff_t>();
    return spsp(g, mask, vs, vt, path, parents);

}

template<class GRAPH, class SUBGRAPH_MASK>
bool 
spsp(
    const GRAPH& g, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const size_t vt,
    std::deque<size_t>& path,
    std::vector<std::ptrdiff_t>& parents// = std::vector<std::ptrdiff_t>()
) {
    path.clear();
    if(!mask.vertex(vs) || !mask.vertex(vt)) {
        return false;
    }
    if(vs == vt) {
        path.push_back(vs);
        return true;
    }
    parents.resize(g.numberOfVertices());
    std::fill(parents.begin(), parents.end(), 0);
    parents[vs] = vs + 1;
    parents[vt] = -static_cast<std::ptrdiff_t>(vt) - 1;
    std::queue<size_t> queues[2];
    queues[0].push(vs);
    queues[1].push(vt);
    for(size_t q = 0; true; q = 1 - q) { // infinite loop, alternating queues
        const size_t numberOfNodesAtFront = queues[q].size();
        for(size_t n = 0; n < numberOfNodesAtFront; ++n) {
            const size_t v = queues[q].front();
            queues[q].pop();
            typename GRAPH::AdjacencyIterator it;
            typename GRAPH::AdjacencyIterator end;
            if(q == 0) {
                it = g.adjacenciesFromVertexBegin(v);
                end = g.adjacenciesFromVertexEnd(v);
            }
            else {
                it = g.adjacenciesToVertexBegin(v);
                end = g.adjacenciesToVertexEnd(v);
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
                        parents[it->vertex()] = -static_cast<std::ptrdiff_t>(v) - 1;
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

/// Search for a shortest path from one to another vertex in a graph with **non-negative edge weights** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param vs Source vertex.
/// \param vt Target vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param path A double-ended queue to which the path is written.
/// \param distance the distance to from the source to the target vertex (if there exists a path).
///     if no path is found, path.size() == 0.
///     the data type of this parameter is used to sum up edge weights.
///
template<
    class GRAPH, 
    class EDGE_WEIGHT_ITERATOR, 
    class T
>
inline void
spsp(
    const GRAPH& g, 
    const size_t vs,
    const size_t vt,
    EDGE_WEIGHT_ITERATOR edgeWeights,
    std::deque<size_t>& path,
    T& distance
) {
    std::vector<T> distances(g.numberOfVertices()); 
    std::vector<size_t> parents(g.numberOfVertices()); 
    spsp(g, DefaultSubgraphMask<>(), vs, vt, edgeWeights, path, distance, 
        distances.begin(), parents.begin()
    );
}

/// Search for a shortest path from one to another vertex in a **subgraph** with **non-negative edge weights** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param vt Target vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param path A double-ended queue to which the path is written.
/// \param distance the distance to from the source to the target vertex (if there exists a path).
///     if no path is found, path.size() == 0.
///     the data type of this parameter is used to sum up edge weights.
///
template<
    class GRAPH, 
    class SUBGRAPH_MASK, 
    class EDGE_WEIGHT_ITERATOR, 
    class T
>
inline void
spsp(
    const GRAPH& g, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const size_t vt,
    EDGE_WEIGHT_ITERATOR edgeWeights,
    std::deque<size_t>& path,
    T& distance
) {
    std::vector<T> distances(g.numberOfVertices()); 
    std::vector<size_t> parents(g.numberOfVertices()); 
    spsp(g, mask, vs, vt, edgeWeights, path, distance,
        distances.begin(), parents.begin()
    );
}

/// Search for shortest paths from a given vertex to every other vertex in an **unweighted** graph using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param vs Source vertex.
/// \param distances Random access iterator pointing to distances
///
template<class GRAPH, class DISTANCE_ITERATOR>
inline void
sssp(
    const GRAPH& g, 
    const size_t vs,
    DISTANCE_ITERATOR distances
) {
    std::vector<size_t> parents(g.numberOfVertices());
    sssp(g, vs, distances, parents.begin());
}

/// Search for shortest paths from a given vertex to every other vertex in an **unweighted** graph using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param vs Source vertex.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
sssp(
    const GRAPH& g, 
    const size_t vs,
    DISTANCE_ITERATOR distances,
    PARENT_ITERATOR parents
) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
    sssp(g, DefaultSubgraphMask<>(), vs, UnitEdgeWeightIterator<Value>(), 
        distances, parents
    );
}

/// Search for shortest paths from a given vertex to every other vertex in an **unweighted** **subgraph** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param distances Random access iterator pointing to distances
///
template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR>
inline void 
sssp(
    const GRAPH& g, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    DISTANCE_ITERATOR distances
) {
    std::vector<size_t> parents(g.numberOfVertices());
    sssp(g, mask, vs, distances, parents.begin());
}

/// Search for shortest paths from a given vertex to every other vertex in an **unweighted** **subgraph** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void 
sssp(
    const GRAPH& g, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    DISTANCE_ITERATOR distances,
    PARENT_ITERATOR parents 
) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
    sssp(g, mask, vs, UnitEdgeWeightIterator<Value>(), distances, parents);
}

/// Search for shortest paths from a given vertex to every other vertex in a graph with **non-negative edge weights** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param vs Source vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<class GRAPH, class EDGE_WEIGHT_ITERATOR, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
sssp(
    const GRAPH& g, 
    const size_t vs,
    const EDGE_WEIGHT_ITERATOR edgeWeights,
    DISTANCE_ITERATOR distances,
    PARENT_ITERATOR parents
) {
    sssp(g, DefaultSubgraphMask<>(), vs, edgeWeights, distances, parents);
}

/// Idle visitor for Dijkstra's algorithm.
///
template<class DISTANCE_ITERATOR, class PARENT_ITERATOR>
struct DijkstraIdleVisitor {
    typedef DISTANCE_ITERATOR Distances;
    typedef PARENT_ITERATOR Parents;

    bool operator()(Distances, Parents, const size_t) const
        { return true; /* continue with the algorithm */ }
};

/// Search for shortest paths from a given vertex to every other vertex in a **subgraph** with **non-negative edge weights** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param distances Random access iterator pointing to distances.
/// \param parents Random access iterator pointing to parent vertices.
///
template<
    class GRAPH, 
    class SUBGRAPH_MASK, 
    class EDGE_WEIGHT_ITERATOR, 
    class DISTANCE_ITERATOR, 
    class PARENT_ITERATOR
>
inline void 
sssp(
    const GRAPH& g, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const EDGE_WEIGHT_ITERATOR edgeWeights,
    DISTANCE_ITERATOR distances,
    PARENT_ITERATOR parents
) {
    typedef DijkstraIdleVisitor<DISTANCE_ITERATOR, PARENT_ITERATOR> Visitor;
    Visitor visitor;
    sssp<GRAPH, SUBGRAPH_MASK, EDGE_WEIGHT_ITERATOR, DISTANCE_ITERATOR, PARENT_ITERATOR, Visitor>(
                g, mask, vs, edgeWeights, distances, parents, visitor
    );
}

/// Search for shortest paths from a given vertex to every other vertex in a **subgraph** with **non-negative edge weights** using Dijkstra's algorithm with a visitor.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param distances Random access iterator pointing to distances.
/// \param parents Random access iterator pointing to parent vertices.
/// \param visitor See DijkstraIdleVisitor.
///
template<
    class GRAPH, 
    class SUBGRAPH_MASK, 
    class EDGE_WEIGHT_ITERATOR, 
    class DISTANCE_ITERATOR, 
    class PARENT_ITERATOR,
    class VISITOR
>
void 
sssp(
    const GRAPH& g, 
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const EDGE_WEIGHT_ITERATOR edgeWeights,
    DISTANCE_ITERATOR distances,
    PARENT_ITERATOR parents,
    VISITOR& visitor 
) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
    typedef typename graph_detail::DijkstraQueueEntry<Value> Entry;

    assert(mask.vertex(vs));  
    const Value infinity = std::numeric_limits<Value>::has_infinity 
        ? std::numeric_limits<Value>::infinity() 
        : std::numeric_limits<Value>::max();
    std::priority_queue<Entry> queue;
    queue.push(vs);
    for(size_t v = 0; v < g.numberOfVertices(); ++v) {
        distances[v] = infinity;
        if(mask.vertex(v)) {
            queue.push(Entry(v, infinity));
        }
    }
    distances[vs] = 0;
    while(!queue.empty()) {
        const size_t v = queue.top().vertex_;
        const bool proceed = visitor(distances, parents, v);
        if(!proceed) {
            return;
        }
        queue.pop();
        if(distances[v] == infinity) {
            return;
        }
        for(typename GRAPH::AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
        it != g.adjacenciesFromVertexEnd(v); ++it) {
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

/// Search for a shortest path from one to another vertex in a **subgraph** with **non-negative edge weights** using Dijkstra's algorithm.
///
/// \param g A graph class such as andres::Graph or andres::Digraph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param vs Source vertex.
/// \param vt Target vertex.
/// \param edgeWeights A random access iterator pointing to positive edge weights.
/// \param path A double-ended queue to which the path is written.
/// \param distance the distance to from the source to the target vertex (if there exists a path).
///     if no path is found, path.size() == 0.
///     the data type of this parameter is used to sum up edge weights.
/// \param distances Random access iterator pointing to distances
/// \param parents Random access iterator pointing to parent vertices
///
template<
    class GRAPH,
    class SUBGRAPH_MASK,
    class EDGE_WEIGHT_ITERATOR,
    class T,
    class DISTANCE_ITERATOR,
    class PARENT_ITERATOR
>
inline void
spsp(
    const GRAPH& g,
    const SUBGRAPH_MASK& mask,
    const size_t vs,
    const size_t vt,
    EDGE_WEIGHT_ITERATOR edgeWeights,
    std::deque<size_t>& path,
    T& distance,
    DISTANCE_ITERATOR distances,
    PARENT_ITERATOR parents
) {
    typedef graph_detail::DijkstraSPSPVisitor<DISTANCE_ITERATOR, PARENT_ITERATOR> Visitor;
    Visitor visitor(vs, vt, path);
    sssp<GRAPH, SUBGRAPH_MASK, EDGE_WEIGHT_ITERATOR, DISTANCE_ITERATOR, PARENT_ITERATOR, Visitor>(
        g, mask, vs, edgeWeights, distances, parents, visitor
    );
    distance = distances[vt];
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_SORTEST_PATHS_HXX
