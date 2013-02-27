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

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

template<class ILP>
void testMulticut() {
    typedef ILP Ilp;
    typedef andres::graph::Graph<> Graph;
    typedef andres::graph::Multicut<Graph, Ilp> Multicut;

    // a simple weighted graph in which an optimal multicut is non-trivial

    Graph graph;
    graph.insertVertices(6);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(0, 3); // 1
    graph.insertEdge(1, 2); // 2
    graph.insertEdge(1, 4); // 3
    graph.insertEdge(2, 5); // 4
    graph.insertEdge(3, 4); // 5
    graph.insertEdge(4, 5); // 6

    std::vector<double> weights(7);
    weights[0] = 5;
    weights[1] = -20;
    weights[2] = 5;
    weights[3] = 5;
    weights[4] = -20;
    weights[5] = 5;
    weights[6] = 5;

    Multicut mc;
    mc.ilp().setVerbosity(false);
    mc.setup(graph, weights);
    mc.solve(10);

    test(mc.label(0) == 0);
    test(mc.label(1) == 1);
    test(mc.label(2) == 0);
    test(mc.label(3) == 1);
    test(mc.label(4) == 1);
    test(mc.label(5) == 0);
    test(mc.label(6) == 0);
}
