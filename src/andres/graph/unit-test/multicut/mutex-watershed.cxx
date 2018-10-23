#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/mutex-watershed.hxx"

inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

void testMutexWatershed() {
    // a simple weighted graph with a unique MWS solution

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
    weights[1] = -21;
    weights[2] = 4;
    weights[3] = 6;
    weights[4] = -20;
    weights[5] = 7;
    weights[6] = 8;

    std::vector<int> edge_labels(graph.numberOfEdges());
    andres::graph::multicut::mutexWatershed(graph, weights, edge_labels);

    test(edge_labels[0] == 1);
    test(edge_labels[1] == 1);
    test(edge_labels[2] == 1);
    test(edge_labels[3] == 0);
    test(edge_labels[4] == 1);
    test(edge_labels[5] == 0);
    test(edge_labels[6] == 0);
}

int main()
{
    testMutexWatershed();
    return 0;
}
