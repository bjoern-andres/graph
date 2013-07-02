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

struct SubgraphMask2 {
    bool vertex(const size_t v) const { return true; }
    bool edge(const size_t e) const { return e != 3; }
};

struct SubgraphMask3 {
    bool vertex(const size_t v) const { return v != 1 && v != 3; }
    bool edge(const size_t e) const { return e != 0 && e != 2 && e!= 4 && e!= 5; }
};

// "incomplete" mask

struct SubgraphMask4 {
    bool vertex(const size_t v) const { return v != 3; }
    bool edge(const size_t e) const { return true; }
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
    
    test(maxFlowPushRelabel.flow(0) == 10);
    test(maxFlowPushRelabel.flow(1) == 5);
    test(maxFlowPushRelabel.flow(2) == 10);
    test(maxFlowPushRelabel.flow(3) == 5);
    test(maxFlowPushRelabel.flow(4) == 5);
    test(maxFlowPushRelabel.flow(5) == 5);
    test(maxFlowPushRelabel.flow(6) == 10);
    test(maxFlowPushRelabel.maxFlow() == 15);
    
    // with subgraph masks
    
    maxFlowPushRelabel(
        digraph,
        SubgraphMask1(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowPushRelabel.flow(0) == 0);
    test(maxFlowPushRelabel.flow(1) == 5);
    test(maxFlowPushRelabel.flow(2) == 0);
    test(maxFlowPushRelabel.flow(3) == 5);
    test(maxFlowPushRelabel.flow(4) == 0);
    test(maxFlowPushRelabel.flow(5) == 0);
    test(maxFlowPushRelabel.flow(6) == 5);
    test(maxFlowPushRelabel.maxFlow() == 5);
    
    maxFlowPushRelabel(
        digraph,
        SubgraphMask2(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowPushRelabel.flow(0) == 10);
    test(maxFlowPushRelabel.flow(1) == 0);
    test(maxFlowPushRelabel.flow(2) == 10);
    test(maxFlowPushRelabel.flow(3) == 0);
    test(maxFlowPushRelabel.flow(4) == 5);
    test(maxFlowPushRelabel.flow(5) == 5);
    test(maxFlowPushRelabel.flow(6) == 5);
    test(maxFlowPushRelabel.maxFlow() == 10);
    
    maxFlowPushRelabel(
        digraph,
        SubgraphMask3(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowPushRelabel.flow(0) == 0);
    test(maxFlowPushRelabel.flow(1) == 5);
    test(maxFlowPushRelabel.flow(2) == 0);
    test(maxFlowPushRelabel.flow(3) == 5);
    test(maxFlowPushRelabel.flow(4) == 0);
    test(maxFlowPushRelabel.flow(5) == 0);
    test(maxFlowPushRelabel.flow(6) == 5);
    test(maxFlowPushRelabel.maxFlow() == 5);
    
    maxFlowPushRelabel.clear();
    
    maxFlowPushRelabel(
        digraph,
        SubgraphMask4(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowPushRelabel.flow(0) == 0);
    test(maxFlowPushRelabel.flow(1) == 5);
    test(maxFlowPushRelabel.flow(2) == 0);
    test(maxFlowPushRelabel.flow(3) == 5);
    test(maxFlowPushRelabel.flow(4) == 0);
    test(maxFlowPushRelabel.flow(5) == 0);
    test(maxFlowPushRelabel.flow(6) == 5);
    test(maxFlowPushRelabel.maxFlow() == 5);
        
    maxFlowPushRelabel.clear();

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
    
    // with subgraph masks
    
    maxFlowEdmondsKarp(
        digraph,
        SubgraphMask1(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowEdmondsKarp.flow(0) == 0);
    test(maxFlowEdmondsKarp.flow(1) == 5);
    test(maxFlowEdmondsKarp.flow(2) == 0);
    test(maxFlowEdmondsKarp.flow(3) == 5);
    test(maxFlowEdmondsKarp.flow(4) == 0);
    test(maxFlowEdmondsKarp.flow(5) == 0);
    test(maxFlowEdmondsKarp.flow(6) == 5);
    test(maxFlowEdmondsKarp.maxFlow() == 5);
    
    maxFlowEdmondsKarp(
        digraph,
        SubgraphMask2(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowEdmondsKarp.flow(0) == 10);
    test(maxFlowEdmondsKarp.flow(1) == 0);
    test(maxFlowEdmondsKarp.flow(2) == 10);
    test(maxFlowEdmondsKarp.flow(3) == 0);
    test(maxFlowEdmondsKarp.flow(4) == 2);
    test(maxFlowEdmondsKarp.flow(5) == 8);
    test(maxFlowEdmondsKarp.flow(6) == 2);
    test(maxFlowEdmondsKarp.maxFlow() == 10);
    
    maxFlowEdmondsKarp(
        digraph,
        SubgraphMask3(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );
    
    test(maxFlowEdmondsKarp.flow(0) == 0);
    test(maxFlowEdmondsKarp.flow(1) == 5);
    test(maxFlowEdmondsKarp.flow(2) == 0);
    test(maxFlowEdmondsKarp.flow(3) == 5);
    test(maxFlowEdmondsKarp.flow(4) == 0);
    test(maxFlowEdmondsKarp.flow(5) == 0);
    test(maxFlowEdmondsKarp.flow(6) == 5);
    test(maxFlowEdmondsKarp.maxFlow() == 5);
    
    maxFlowEdmondsKarp.clear();
    
    maxFlowEdmondsKarp(
        digraph,
        SubgraphMask4(),
        edgeWeights.begin(),
        sourceVertexIndex,
        sinkVertexIndex
    );

    test(maxFlowEdmondsKarp.flow(1) == 5);
    test(maxFlowEdmondsKarp.flow(2) == 0);
    test(maxFlowEdmondsKarp.flow(3) == 5);
    test(maxFlowEdmondsKarp.flow(4) == 0);
    test(maxFlowEdmondsKarp.flow(5) == 0);
    test(maxFlowEdmondsKarp.flow(6) == 5);
    test(maxFlowEdmondsKarp.maxFlow() == 5);

    maxFlowEdmondsKarp.clear();

}

int main() {
    testPushRelabel();
    testEdmondsKarp();

    return 0;
}
