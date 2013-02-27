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
#include "andres/graph/shortest-paths.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

int main() {
    // spsp, undirected graph
    {
        andres::graph::Graph<> g(9);
        g.insertEdge(0, 1);
        g.insertEdge(1, 2);
        g.insertEdge(2, 3);
        g.insertEdge(0, 4);
        g.insertEdge(3, 4);
        g.insertEdge(0, 5);
        g.insertEdge(5, 3);
        g.insertEdge(0, 6);
        g.insertEdge(6, 7);
        g.insertEdge(7, 3);

        std::deque<size_t> path;

        bool found = andres::graph::spsp(g, 0, 3, path);
        test(found == true);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 4);
        test(path[2] == 3);

        found = andres::graph::spsp(g, 0, 8, path);
        test(found == false);
        test(path.size() == 0);

        found = andres::graph::spsp(g, 0, 0, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 0);

        found = andres::graph::spsp(g, 8, 8, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 8);
    }

    // spsp, directed graph
    {
        andres::graph::Digraph<> g(9);
        g.insertEdge(0, 1); // 0
        g.insertEdge(1, 2); // 1
        g.insertEdge(2, 3); // 2
        g.insertEdge(0, 4); // 3
        g.insertEdge(3, 4); // 4
        g.insertEdge(0, 5); // 5
        g.insertEdge(5, 3); // 6
        g.insertEdge(0, 6); // 7
        g.insertEdge(6, 7); // 8
        g.insertEdge(7, 3); // 9

        std::deque<size_t> path;

        bool found = andres::graph::spsp(g, 0, 3, path);   
        test(found == true);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);

        found = andres::graph::spsp(g, 3, 0, path);
        test(found == false);
        test(path.size() == 0);

        found = andres::graph::spsp(g, 0, 8, path);
        test(found == false);
        test(path.size() == 0);

        found = andres::graph::spsp(g, 0, 0, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 0);

        found = andres::graph::spsp(g, 8, 8, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 8);

        // with subgraph mask
        struct SubgraphMask {
            bool vertex(const size_t v) const { return v != 1; }
            bool edge(const size_t e) const { return e < 3 || e > 6; }
        };

        found = andres::graph::spsp(g, SubgraphMask(), 0, 3, path);
        test(found == true);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 6);
        test(path[2] == 7);
        test(path[3] == 3);
    }

    return 0;
}