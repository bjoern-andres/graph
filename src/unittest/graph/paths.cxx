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
#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/paths.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

struct SubgraphMask {
    bool vertex(const size_t v) const { return true; }
    bool edge(const size_t e) const { return e != 7; }
};

int main() {
    typedef andres::graph::Graph<> Graph;

    Graph graph(7);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(1, 2); // 1
    graph.insertEdge(2, 3); // 2
    graph.insertEdge(3, 4); // 3
    graph.insertEdge(4, 5); // 4
    graph.insertEdge(1, 6); // 5
    graph.insertEdge(6, 2); // 6
    graph.insertEdge(3, 5); // 7

    {
        size_t path[] =  {0, 1, 2, 3, 4, 5};

        std::pair<bool, size_t> p = andres::graph::findChord(graph, path, path + 6);
        test(p.first == true);
        test(p.second == 7);

        p = andres::graph::findChord(graph, path, path + 5);
        test(p.first == false);

        p = andres::graph::findChord(graph, path + 3, path + 6);
        test(p.first == true);
        test(p.second == 7);

        p = andres::graph::findChord(graph, SubgraphMask(), path + 3, path + 6, false);
        test(p.first == false);
    }
    {
        size_t path[] = {1, 6, 2};

        std::pair<bool, size_t> p = andres::graph::findChord(graph, path, path + 3);
        test(p.first == true);
        test(p.second == 1);

        p = andres::graph::findChord(graph, path, path + 3, true);
        test(p.first == false);
    }

    return 0;
}
