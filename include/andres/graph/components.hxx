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
#ifndef ANDRES_GRAPH_COMPONENTS_HXX
#define ANDRES_GRAPH_COMPONENTS_HXX

#include <vector>
#include <queue>
#include <algorithm> // std::fill

#include "andres/partition.hxx"

namespace andres {
namespace graph {

/// Connected component labeling by breadth-first-search (labels start at 0).
template<class GRAPH>
struct ComponentsBySearch {
    typedef GRAPH Graph;
    
    ComponentsBySearch();
    size_t build(const Graph&);
    template<class SUBGRAPH_MASK>
        size_t build(const Graph&, const SUBGRAPH_MASK&);
    bool areConnected(const size_t, const size_t) const;

    std::vector<size_t> labels_;
};

/// Connected component labeling using disjoint sets (labels start at 0).
template<class GRAPH>
struct ComponentsByPartition {
    typedef GRAPH Graph;
    
    ComponentsByPartition();
    size_t build(const Graph&);
    template<class SUBGRAPH_MASK>
        size_t build(const Graph&, const SUBGRAPH_MASK&);
    bool areConnected(const size_t, const size_t) const;

    andres::Partition<size_t> partition_;
};

/// Connected component labeling by breadth-first-search (labels start at 0).
///
/// \param graph Graph.
/// \param labeling Random access iterator to a container that has as many
///        entries as there are vertices in the graph; all entries need to
///        be initialized as 0
///
template<class GRAPH, class ITERATOR>
inline size_t
labelComponents(
    const GRAPH& graph,
    ITERATOR labeling
) {
    return labelComponents(graph, DefaultSubgraphMask<>(), labeling);
}

/// Connected component labeling by breadth-first-search (labels start at 0).
///
/// \param graph Graph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param labeling Random access iterator to a container that has as many
///        entries as there are vertices in the graph; all entries need to
///        be initialized as 0
///
template<class GRAPH, class SUBGRAPH_MASK, class ITERATOR>
size_t
labelComponents(
    const GRAPH& graph,
    const SUBGRAPH_MASK& mask,
    ITERATOR labeling
) {
    size_t label = 0;
    std::vector<bool> visited(graph.numberOfVertices(), false);
    std::queue<size_t> queue;
    for(size_t v = 0; v < graph.numberOfVertices(); ++v) {
        if(mask.vertex(v)) {
            if(!visited[v]) {
                labeling[v] = label; // label
                queue.push(v);
                visited[v] = true;
                while(!queue.empty()) {
                    size_t w = queue.front();
                    queue.pop();
                    for(GRAPH::AdjacencyIterator it = graph.adjacenciesFromVertexBegin(w);
                    it != graph.adjacenciesFromVertexEnd(w); ++it) {
                        if(mask.edge(it->edge()) 
                        && mask.vertex(it->vertex()) 
                        && !visited[it->vertex()]) {
                            labeling[it->vertex()] = label; // label
                            queue.push(it->vertex());
                            visited[it->vertex()] = true;
                        }
                    }
                }
                label++;
            }
        }
        else {
            labeling[v] = 0;
        }
    }
    return label;
}

template<class GRAPH>
inline 
ComponentsBySearch<GRAPH>::ComponentsBySearch()
:   labels_()
{}

template<class GRAPH>
inline size_t
ComponentsBySearch<GRAPH>::build(
    const Graph& graph
) {
    return build(graph, DefaultSubgraphMask<>());
}

template<class GRAPH>
template<class SUBGRAPH_MASK>
inline size_t 
ComponentsBySearch<GRAPH>::build(
    const Graph& graph,
    const SUBGRAPH_MASK& mask
) {
    labels_.resize(graph.numberOfVertices());
    return labelComponents(graph, mask, labels_.begin());
}

template<class GRAPH>
inline bool 
ComponentsBySearch<GRAPH>::areConnected(
    const size_t vertex0, 
    const size_t vertex1
) const {
    return labels_[vertex0] == labels_[vertex1];
}

template<class GRAPH>
inline 
ComponentsByPartition<GRAPH>::ComponentsByPartition()
:   partition_()
{}

template<class GRAPH>
inline size_t
ComponentsByPartition<GRAPH>::build(
    const Graph& graph
) {
    return build(graph, DefaultSubgraphMask<>());
}

template<class GRAPH>
template<class SUBGRAPH_MASK>
inline size_t 
ComponentsByPartition<GRAPH>::build(
    const Graph& graph,
    const SUBGRAPH_MASK& mask
) {
    partition_.assign(graph.numberOfVertices());
    for(size_t edge = 0; edge < graph.numberOfEdges(); ++edge) {
        if(mask.edge(edge)) {
            const size_t v0 = graph.vertexOfEdge(edge, 0);
            const size_t v1 = graph.vertexOfEdge(edge, 1);
            if(mask.vertex(v0) && mask.vertex(v1)) {
                partition_.merge(v0, v1);
            }
        }
    }
    return partition_.numberOfSets();
}

template<class GRAPH>
inline bool 
ComponentsByPartition<GRAPH>::areConnected(
    const size_t vertex0, 
    const size_t vertex1
) const {
    return partition_.find(vertex0) == partition_.find(vertex1);
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_COMPONENTS_HXX
