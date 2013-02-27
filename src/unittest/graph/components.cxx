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
#include "andres/graph/components.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

template<class COMPONENTS>
void testComponents() {
    typedef COMPONENTS Components;
    typedef Components::Graph Graph;

    Graph graph(10);
    graph.insertEdge(0, 1);
    graph.insertEdge(2, 3);
    graph.insertEdge(3, 4);
    graph.insertEdge(5, 6);
    graph.insertEdge(5, 7);
    graph.insertEdge(6, 8);
    graph.insertEdge(7, 8);

    Components components;
    components.build(graph);

    test(components.areConnected(0, 1));
    test(components.areConnected(2, 3));
    test(components.areConnected(3, 4));
    test(components.areConnected(5, 6));
    test(components.areConnected(6, 7));
    test(components.areConnected(7, 8));
    test(!components.areConnected(0, 2));
    test(!components.areConnected(0, 5));
    test(!components.areConnected(0, 9));
    test(!components.areConnected(2, 5));
    test(!components.areConnected(2, 9));
    test(!components.areConnected(5, 9));

    struct SubgraphMask {
        bool vertex(const size_t v) const { return v != 3; }
        bool edge(const size_t e) const { return e != 0; }
    };

    components.build(graph, SubgraphMask());

    test(components.areConnected(5, 6));
    test(components.areConnected(6, 7));
    test(components.areConnected(7, 8));
    test(!components.areConnected(0, 1));
    test(!components.areConnected(0, 2));
    test(!components.areConnected(0, 4));
    test(!components.areConnected(0, 5));
    test(!components.areConnected(0, 9));
    test(!components.areConnected(1, 2));
    test(!components.areConnected(1, 4));
    test(!components.areConnected(1, 5));
    test(!components.areConnected(1, 9));
    test(!components.areConnected(2, 4));
    test(!components.areConnected(2, 5));
    test(!components.areConnected(2, 9));
    test(!components.areConnected(4, 5));
    test(!components.areConnected(4, 9));
    test(!components.areConnected(5, 9));
}

int main() {
    typedef andres::graph::Graph<> Graph;
    typedef andres::graph::ComponentsBySearch<Graph> ComponentsBySearch;
    typedef andres::graph::ComponentsByPartition<Graph> ComponentsByPartition;

    testComponents<ComponentsBySearch>();
    testComponents<ComponentsByPartition>();

    return 0;
}
