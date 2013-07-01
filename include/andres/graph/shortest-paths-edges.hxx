//
//  shortest-paths-edges.hxx
//  graph
//
//  Created by Duligur Ibeling on 6/18/13.
//
//
#pragma once
#ifndef ANDRES_GRAPH_SORTEST_PATHS_EDGES_HXX
#define ANDRES_GRAPH_SORTEST_PATHS_EDGES_HXX

#include <limits>
#include <deque>
#include <queue>
#include <vector>

#include "andres/graph/graph.hxx"
#include "andres/graph/shortest-paths.hxx"

namespace andres {
namespace graph {

template<class GRAPH>
inline bool
spspEdges(
    const GRAPH&,
    const size_t,
    const size_t,
    std::deque<size_t>&,
    std::vector<std::ptrdiff_t>&
);
    
template<class GRAPH>
inline bool
spspEdges(
    const GRAPH&,
    const size_t,
    const size_t,
    std::deque<size_t>&
);
    
template<class GRAPH, class SUBGRAPH_MASK>
inline bool
spspEdges(
    const GRAPH&,
    const SUBGRAPH_MASK&,
    const size_t,
    const size_t,
    std::deque<size_t>&
);
    
template<class GRAPH, class SUBGRAPH_MASK>
bool
spspEdges(
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
spspEdges(
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
spspEdges(
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
ssspEdges(
     const GRAPH&,
     const size_t,
     DISTANCE_ITERATOR
);

template<class GRAPH, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const size_t,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR
);

template<class GRAPH, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const size_t,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR,
	 PARENT_ITERATOR
);

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const SUBGRAPH_MASK&,
     const size_t,
     DISTANCE_ITERATOR
);

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const SUBGRAPH_MASK&,
     const size_t,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR
);

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const SUBGRAPH_MASK&,
     const size_t,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR,
	 PARENT_ITERATOR
);

template<class GRAPH, class EDGE_WEIGHT_ITERATOR, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const size_t,
     const EDGE_WEIGHT_ITERATOR,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR
);

template<class GRAPH, class EDGE_WEIGHT_ITERATOR, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH&,
     const size_t,
     const EDGE_WEIGHT_ITERATOR,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR,
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
ssspEdges(
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
class PARENT_ITERATOR
>
inline void
ssspEdges(
     const GRAPH&,
     const SUBGRAPH_MASK&,
     const size_t,
     const EDGE_WEIGHT_ITERATOR,
     DISTANCE_ITERATOR,
     PARENT_ITERATOR,
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
ssspEdges(
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
spspEdges(
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

template<class GRAPH>
inline bool
spspEdges(
          const GRAPH& g,
          const size_t vs,
          const size_t vt,
          std::deque<size_t>& path,
          std::vector<std::ptrdiff_t>& parents
          ) {
    return spspEdges(g, DefaultSubgraphMask<>(), vs, vt, path, parents);
}
    
template<class GRAPH>
inline bool
spspEdges(
          const GRAPH& g,
          const size_t vs,
          const size_t vt,
          std::deque<size_t>& path
          ) {
    std::vector<std::ptrdiff_t> parents = std::vector<std::ptrdiff_t>();
    return spspEdges(g, DefaultSubgraphMask<>(), vs, vt, path, parents);
}

template<class GRAPH, class SUBGRAPH_MASK>
inline bool
spspEdges(
          const GRAPH& g,
          const SUBGRAPH_MASK& mask,
          const size_t vs,
          const size_t vt,
          std::deque<size_t>& path
) {
    std::vector<std::ptrdiff_t> parents = std::vector<std::ptrdiff_t>();
    return spspEdges(g, mask, vs, vt, path, parents);
}
    
template<class GRAPH, class SUBGRAPH_MASK>
bool
spspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     const size_t vt,
     std::deque<size_t>& path, // sequence of edges
     std::vector<std::ptrdiff_t>& parents // sequence of edges
) {
    path.clear();
    if(!mask.vertex(vs) || !mask.vertex(vt)) {
        return false;
    }
    if(vs == vt) {
        return true;
    }
    for (typename GRAPH::AdjacencyIterator i = g.adjacenciesFromVertexBegin(vs); i < g.adjacenciesFromVertexEnd(vs) ; ++i) {
        if (i->vertex() == vt && mask.edge(i->edge())) {
            path.push_front(i->edge());
            return true;
        }
    }
    parents.resize(g.numberOfEdges());
    std::fill(parents.begin(), parents.end(), 0);
    std::queue<size_t> queues[2];
    for (typename GRAPH::AdjacencyIterator i = g.adjacenciesFromVertexBegin(vs); i < g.adjacenciesFromVertexEnd(vs) ; ++i) {
        if (mask.edge(i->edge())) {
            queues[0].push(i->edge());
            parents[i->edge()] = i->edge() + 1;
        }        
    }
    for (typename GRAPH::AdjacencyIterator i = g.adjacenciesToVertexBegin(vt); i < g.adjacenciesToVertexEnd(vt) ; ++i) {
        if (mask.edge(i->edge())) {
        queues[1].push(i->edge());
        parents[i->edge()] = -static_cast<std::ptrdiff_t>(i->edge()) - 1;
        }
    }
    for(size_t q = 0; true; q = 1 - q) { // infinite loop, alternating queues
        const size_t numberOfEdgesAtFront = queues[q].size();
        for(size_t n = 0; n < numberOfEdgesAtFront; ++n) {
            const size_t e = queues[q].front();
            queues[q].pop();
            typename GRAPH::AdjacencyIterator it;
            typename GRAPH::AdjacencyIterator end;
            if(q == 0) {
                it = g.adjacenciesFromVertexBegin(g.vertexOfEdge(e, 1));
                end = g.adjacenciesFromVertexEnd(g.vertexOfEdge(e, 1));
            }
            else {
                it = g.adjacenciesToVertexBegin(g.vertexOfEdge(e, 0));
                end = g.adjacenciesToVertexEnd(g.vertexOfEdge(e, 0));
            }
            for(; it != end; ++it) {
                if(!mask.edge(it->edge()) || !mask.vertex(it->vertex())) {
                    continue;
                }
                if(parents[it->edge()] < 0 && q == 0) {
                    graph_detail::spspHelper(parents, e, it->edge(), path);
                    return true;
                }
                else if(parents[it->edge()] > 0 && q == 1) {
                    graph_detail::spspHelper(parents, it->edge(), e, path);
                    return true;
                }
                else if(parents[it->edge()] == 0) {
                    if(q == 0) {
                        parents[it->edge()] = e + 1;
                    }
                    else {
                        parents[it->edge()] = -static_cast<std::ptrdiff_t>(e) - 1;
                    }
                    queues[q].push(it->edge());
                }
            }
        }
        if(queues[0].empty() && queues[1].empty()) {
            return false;
        }
    }

}
    
template<
class GRAPH,
class EDGE_WEIGHT_ITERATOR,
class T
>
inline void
spspEdges(
     const GRAPH& g,
     const size_t vs,
     const size_t vt,
     EDGE_WEIGHT_ITERATOR edgeWeights,
     std::deque<size_t>& path,
     T& distance
) {
    std::vector<T> distances(g.numberOfVertices());
    std::vector<size_t> parents(g.numberOfVertices());
	std::vector<size_t> parentsEdges(g.numberOfVertices());
    spspEdges(g, DefaultSubgraphMask<>(), vs, vt, edgeWeights, path, distance,
         distances.begin(), parents.begin(), parentsEdges.begin());
}
    
template<
class GRAPH,
class SUBGRAPH_MASK,
class EDGE_WEIGHT_ITERATOR,
class T
>
inline void
spspEdges(
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
	std::vector<size_t> parentsEdges(g.numberOfVertices());
    spspEdges(g, mask, vs, vt, edgeWeights, path, distance,
         distances.begin(), parents.begin(), parentsEdges.begin());
}
 
template<class GRAPH, class DISTANCE_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const size_t vs,
     DISTANCE_ITERATOR distances
) {
    std::vector<size_t> parents(g.numberOfVertices());
    ssspEdges(g, vs, distances, parents.begin());
}
    
template<class GRAPH, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const size_t vs,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents
) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
	std::vector<size_t> parentsEdges(g.numberOfVertices());
    ssspEdges(g, DefaultSubgraphMask<>(), vs, UnitEdgeWeightIterator<Value>(),
         distances, parents, parentsEdges.begin());
}

template<class GRAPH, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const size_t vs,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents,
	 PARENT_ITERATOR parentsEdges
) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
    ssspEdges(g, DefaultSubgraphMask<>(), vs, UnitEdgeWeightIterator<Value>(),
         distances, parents, parentsEdges);
}

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     DISTANCE_ITERATOR distances
) {
    std::vector<size_t> parents(g.numberOfVertices());
    ssspEdges(g, mask, vs, distances, parents.begin());
}

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents
     ) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
	std::vector<size_t> parentsEdges(g.numberOfVertices());
    ssspEdges(g, mask, vs, UnitEdgeWeightIterator<Value>(), distances, parentsEdges.begin());
}

template<class GRAPH, class SUBGRAPH_MASK, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents,
	 PARENT_ITERATOR parentsEdges
) {
    typedef typename std::iterator_traits<DISTANCE_ITERATOR>::value_type Value;
    ssspEdges(g, mask, vs, UnitEdgeWeightIterator<Value>(), distances, parents, parentsEdges);
}
    
template<class GRAPH, class EDGE_WEIGHT_ITERATOR, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const size_t vs,
     const EDGE_WEIGHT_ITERATOR edgeWeights,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents
     ) {
    ssspEdges(g, DefaultSubgraphMask<>(), vs, edgeWeights, distances, parents);
}

template<class GRAPH, class EDGE_WEIGHT_ITERATOR, class DISTANCE_ITERATOR, class PARENT_ITERATOR>
inline void
ssspEdges(
     const GRAPH& g,
     const size_t vs,
     const EDGE_WEIGHT_ITERATOR edgeWeights,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents,
     PARENT_ITERATOR parentsEdges
     ) {
    ssspEdges(g, DefaultSubgraphMask<>(), vs, edgeWeights, distances, parents, parentsEdges);
}

template<
class GRAPH,
class SUBGRAPH_MASK,
class EDGE_WEIGHT_ITERATOR,
class DISTANCE_ITERATOR,
class PARENT_ITERATOR
>
inline void
ssspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     const EDGE_WEIGHT_ITERATOR edgeWeights,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents
     ) {
    typedef DijkstraIdleVisitor<DISTANCE_ITERATOR, PARENT_ITERATOR> Visitor;
    Visitor visitor;
	std::vector<size_t> parentsEdges(g.numberOfVertices());
    ssspEdges<GRAPH, SUBGRAPH_MASK, EDGE_WEIGHT_ITERATOR, DISTANCE_ITERATOR, PARENT_ITERATOR, Visitor>(g, mask, vs, edgeWeights, distances, parents, parentsEdges.begin(), visitor);
}

template<
class GRAPH,
class SUBGRAPH_MASK,
class EDGE_WEIGHT_ITERATOR,
class DISTANCE_ITERATOR,
class PARENT_ITERATOR
>
inline void
ssspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     const EDGE_WEIGHT_ITERATOR edgeWeights,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents,
	 PARENT_ITERATOR parentsEdges
     ) {
    typedef DijkstraIdleVisitor<DISTANCE_ITERATOR, PARENT_ITERATOR> Visitor;
    Visitor visitor;
    ssspEdges<GRAPH, SUBGRAPH_MASK, EDGE_WEIGHT_ITERATOR, DISTANCE_ITERATOR, PARENT_ITERATOR, Visitor>(
		g, mask, vs, edgeWeights, distances, parents, parentsEdges, visitor);
}

template<
class GRAPH,
class SUBGRAPH_MASK,
class EDGE_WEIGHT_ITERATOR,
class DISTANCE_ITERATOR,
class PARENT_ITERATOR,
class VISITOR
>
void
ssspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     const EDGE_WEIGHT_ITERATOR edgeWeights,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents,
     PARENT_ITERATOR parentsEdges,
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
        const bool proceed = visitor(distances, parents, parentsEdges, v);
        if(!proceed) {
            // return;
			break;
        }
        queue.pop();
        if(distances[v] == infinity) {
            // return;
			break;
        }
        for(typename GRAPH::AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
            it != g.adjacenciesFromVertexEnd(v); ++it) {
            if(mask.vertex(it->vertex()) && mask.edge(it->edge())) {
                const Value alternativeDistance = distances[v] + edgeWeights[it->edge()];
                if(alternativeDistance < distances[it->vertex()]) {
                    distances[it->vertex()] = alternativeDistance;
                    parentsEdges[it->vertex()] = it->edge();
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

template<
class GRAPH,
class SUBGRAPH_MASK,
class EDGE_WEIGHT_ITERATOR,
class T,
class DISTANCE_ITERATOR,
class PARENT_ITERATOR
>
inline void
spspEdges(
     const GRAPH& g,
     const SUBGRAPH_MASK& mask,
     const size_t vs,
     const size_t vt,
     EDGE_WEIGHT_ITERATOR edgeWeights,
     std::deque<size_t>& path,
     T& distance,
     DISTANCE_ITERATOR distances,
     PARENT_ITERATOR parents,
	 PARENT_ITERATOR parentsEdges
)
{
    typedef graph_detail::DijkstraSPSPVisitor<DISTANCE_ITERATOR, PARENT_ITERATOR> Visitor;
    Visitor visitor(vs, vt, path);
    ssspEdges<GRAPH, SUBGRAPH_MASK, EDGE_WEIGHT_ITERATOR, DISTANCE_ITERATOR, PARENT_ITERATOR, Visitor>(g, mask, vs, edgeWeights, distances, parents, parentsEdges, visitor);
    distance = distances[vt];
}
    
} // namespace graph
} // namespace andres

#endif
