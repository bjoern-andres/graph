/// Copyright (c) 2013 by Mark Matten
/// 
/// This software was developed by Mark Matten.
/// Enquiries shall be directed to markmatten@gmail.com
///
/// Redistribution and use in source and binary forms, with or without 
/// modification, are permitted provided that the following conditions are met:
///
/// - Redistributions of source code must retain the above copyright notice,
///   this list of conditions and the following disclaimer.
/// - Redistributions in binary form must reproduce the above copyright notice, 
///   this list of conditions and the following disclaimer in the documentation
///   and/or other materials provided with the distribution.
/// - The name of the author must not be used to endorse or promote products 
///   derived from this software without specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
/// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
/// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
/// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
/// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
/// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
/// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
/// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
/// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///
#include <stdexcept>
#include <vector>

#include "andres/graph/graph.hxx"
#include "andres/graph/max-flow.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

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

	test(maxFlowPushRelabel.flow(0) == 10);
	test(maxFlowPushRelabel.flow(1) == 5);
	test(maxFlowPushRelabel.flow(2) == 10);
	test(maxFlowPushRelabel.flow(3) == 5);
	test(maxFlowPushRelabel.flow(4) == 5);
	test(maxFlowPushRelabel.flow(5) == 5);
	test(maxFlowPushRelabel.flow(6) == 10);
	test(maxFlowPushRelabel.maxFlow() == 15);

	maxFlowPushRelabel.clear(); // just to see if it works
}

int main() {
	testPushRelabel();

	return 0;
}
