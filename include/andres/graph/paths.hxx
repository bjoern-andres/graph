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
#ifndef ANDRES_GRAPH_PATHS_HXX
#define ANDRES_GRAPH_PATHS_HXX

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
inline std::pair<bool, size_t> 
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
inline std::pair<bool, size_t> 
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
        std::pair<bool, size_t> p = graph.findEdge(*it, *it2);
        if(p.first && mask.edge(p.second)) {
            return p;
        }
    }
    return std::pair<bool, size_t>(false, 0);
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_PATHS_HXX
