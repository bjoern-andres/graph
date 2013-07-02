// Copyright (c) 2013 by Mark Matten, Duligur Ibeling.
// 
// This software was developed by Mark Matten.
// Enquiries shall be directed to markmatten@gmail.com
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
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
#include <vector>

#include "andres/graph/graph.hxx"
#include "andres/graph/max-flow.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

struct SubgraphMask1 {
    bool vertex(const size_t v) const { return true; }
    bool edge(const size_t e) const { return e != 2; }
};

void testPushRelabel() {
	typedef andres::graph::Digraph<> Digraph;
	typedef double Flow;
	typedef andres::graph::MaxFlowPushRelabel<Digraph, Flow> MaxFlowPushRelabel;
	
	// define graph
	Digraph digraph;
	digraph.insertVertices(6);
	digraph.insertEdge(0, 1);
	digraph.insertEdge(0, 2);
	digraph.insertEdge(1, 3);
	digraph.insertEdge(2, 4);
	digraph.insertEdge(3, 4);
	digraph.insertEdge(3, 5);
	digraph.insertEdge(4, 5);

	// define edge weights
	std::vector<double> edgeWeights(digraph.numberOfEdges());
	edgeWeights[0] = 20;
	edgeWeights[1] = 15;
	edgeWeights[2] = 10;
	edgeWeights[3] = 5;
	edgeWeights[4] = 5;
	edgeWeights[5] = 8;
	edgeWeights[6] = 10;

	// define source and sink vertex
	size_t sourceVertexIndex = 0;
	size_t sinkVertexIndex = 5;

	MaxFlowPushRelabel maxFlowPushRelabel(
		digraph, 
		edgeWeights.begin(), 
		sourceVertexIndex, 
		sinkVertexIndex
	);
    
    MaxFlowPushRelabel maxFlowPushRelabel1(digraph, SubgraphMask1(), edgeWeights.begin(), sourceVertexIndex, sinkVertexIndex);
        
	maxFlowPushRelabel.clear(); // just to see if it works

}

void testEdmondsKarp() {
	typedef andres::graph::Digraph<> Digraph;
	typedef double Flow;
	typedef andres::graph::MaxFlowEdmondsKarp<Digraph, Flow> MaxFlowEdmondsKarp;
	
	// define graph
	Digraph digraph;
	digraph.insertVertices(6);
	digraph.insertEdge(0, 1);
	digraph.insertEdge(0, 2);
	digraph.insertEdge(1, 3);
	digraph.insertEdge(2, 4);
	digraph.insertEdge(3, 4);
	digraph.insertEdge(3, 5);
	digraph.insertEdge(4, 5);
    
	// define edge weights
	std::vector<double> edgeWeights(digraph.numberOfEdges());
	edgeWeights[0] = 20;
	edgeWeights[1] = 15;
	edgeWeights[2] = 10;
	edgeWeights[3] = 5;
	edgeWeights[4] = 5;
	edgeWeights[5] = 8;
	edgeWeights[6] = 10;
    
	// define source and sink vertex
	size_t sourceVertexIndex = 0;
	size_t sinkVertexIndex = 5;
    
	MaxFlowEdmondsKarp maxFlowEdmondsKarp(
        digraph,
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowEdmondsKarp.flow(0) == 10);
	test(maxFlowEdmondsKarp.flow(1) == 5);
	test(maxFlowEdmondsKarp.flow(2) == 10);
	test(maxFlowEdmondsKarp.flow(3) == 5);
	test(maxFlowEdmondsKarp.flow(4) == 2);
	test(maxFlowEdmondsKarp.flow(5) == 8);
	test(maxFlowEdmondsKarp.flow(6) == 7);
	test(maxFlowEdmondsKarp.maxFlow() == 15);
    
    // with subgraph mask
    
    MaxFlowEdmondsKarp maxFlowEdmondsKarp1(
        digraph,
        SubgraphMask1(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowEdmondsKarp1.flow(0) == 0);
	test(maxFlowEdmondsKarp1.flow(1) == 5);
	test(maxFlowEdmondsKarp1.flow(2) == 0);
	test(maxFlowEdmondsKarp1.flow(3) == 5);
	test(maxFlowEdmondsKarp1.flow(4) == 0);
	test(maxFlowEdmondsKarp1.flow(5) == 0);
	test(maxFlowEdmondsKarp1.flow(6) == 5);
    test(maxFlowEdmondsKarp1.maxFlow() == 5);
    
    maxFlowEdmondsKarp.clear();
    maxFlowEdmondsKarp1.clear();
    
}

int main() {
	testPushRelabel();
    // testEdmondsKarp();

	return 0;
}
