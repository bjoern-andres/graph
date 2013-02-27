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

#include <deque>
#include <queue>
#include <vector>

#include "andres/graph/graph.hxx" // DefaultSubgraphMask

namespace andres {
namespace graph {

// \cond SUPPRESS_DOXYGEN
namespace detail {

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

} // namespace detail
// \endcond

/// Search for a shortest path connecting a single pair of vertices in an unweighted graph.
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

/// Search for a shortest path connecting a single pair of vertices in an unweighted SUBgraph.
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
                    detail::spspHelper(parents, v, it->vertex(), path);
                    assert(path[0] == vs);
                    assert(path.back() == vt);
                    return true;
                }
                else if(parents[it->vertex()] > 0 && q == 1) {
                    detail::spspHelper(parents, it->vertex(), v, path);
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

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_SORTEST_PATHS_HXX
