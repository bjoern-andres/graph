#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

void testMulticut() {
    // a simple weighted graph in which an optimal multicut is non-trivial

    andres::graph::Graph<> graph;
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

    std::vector<char> edge_labels(graph.numberOfEdges());
    andres::graph::multicut::greedyAdditiveEdgeContraction(graph, weights, edge_labels);

    test(edge_labels[0] == 0);
    test(edge_labels[1] == 1);
    test(edge_labels[2] == 0);
    test(edge_labels[3] == 1);
    test(edge_labels[4] == 1);
    test(edge_labels[5] == 0);
    test(edge_labels[6] == 0);
}

int main()
{
    testMulticut();

    return 0;
}
