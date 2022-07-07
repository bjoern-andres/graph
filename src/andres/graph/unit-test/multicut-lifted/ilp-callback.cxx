#include <stdexcept>
#include <vector>
#include <iostream>
#include <sstream>

#include <andres/ilp/gurobi-callback.hxx>
#include <andres/graph/graph.hxx>
#include <andres/graph/multicut-lifted/ilp-callback.hxx>

inline void test(bool pred) {
    if(!pred) throw std::logic_error("test failed.");
}

void test_3() {
    andres::graph::Graph<> original_graph(3);
    original_graph.insertEdge(0, 1);
    original_graph.insertEdge(1, 2);

    andres::graph::Graph<> lifted_graph(original_graph.numberOfVertices());
    lifted_graph.insertEdge(0, 1);
    lifted_graph.insertEdge(1, 2);
    lifted_graph.insertEdge(0, 2);

    std::vector<double> edge_values(lifted_graph.numberOfEdges());
    edge_values[0] = -10;
    edge_values[1] = -10;
    edge_values[2] = 10;

    std::vector<char> edge_labels(lifted_graph.numberOfEdges());
    andres::graph::multicut_lifted::ilp_callback<andres::ilp::GurobiCallback>(original_graph, lifted_graph, edge_values, edge_labels, edge_labels);

    test(edge_labels[0] == 1);
    test(edge_labels[1] == 1);
    test(edge_labels[2] == 1);
}

void test_4() {
    andres::graph::Graph<> original_graph(4);
    original_graph.insertEdge(0, 1);
    original_graph.insertEdge(0, 3);
    original_graph.insertEdge(1, 2);
    original_graph.insertEdge(2, 3);

    andres::graph::Graph<> lifted_graph(original_graph.numberOfVertices());
    lifted_graph.insertEdge(0, 1);
    lifted_graph.insertEdge(0, 3);
    lifted_graph.insertEdge(1, 2);
    lifted_graph.insertEdge(2, 3);
    lifted_graph.insertEdge(0, 2);

    std::vector<double> edge_values(lifted_graph.numberOfEdges());
    edge_values[0] = -10;
    edge_values[1] = -2;
    edge_values[2] = -10;
    edge_values[3] = -4;
    edge_values[4] = 1;

    std::vector<char> edge_labels(lifted_graph.numberOfEdges());
    andres::graph::multicut_lifted::ilp_callback<andres::ilp::GurobiCallback>(original_graph, lifted_graph, edge_values, edge_labels, edge_labels);

    test(edge_labels[0] == 1);
    test(edge_labels[1] == 1);
    test(edge_labels[2] == 1);
    test(edge_labels[3] == 1);
    test(edge_labels[4] == 1);
}

void test_8() {
    andres::graph::Graph<> original_graph(8);
    original_graph.insertEdge(0, 1);
    original_graph.insertEdge(1, 2);
    original_graph.insertEdge(2, 3);
    original_graph.insertEdge(4, 5);
    original_graph.insertEdge(5, 6);
    original_graph.insertEdge(6, 7);
    original_graph.insertEdge(0, 4);
    original_graph.insertEdge(1, 5);
    original_graph.insertEdge(2, 6);
    original_graph.insertEdge(3, 7);

    andres::graph::Graph<> lifted_graph(original_graph.numberOfVertices());
    lifted_graph.insertEdge(0, 1);
    lifted_graph.insertEdge(1, 2);
    lifted_graph.insertEdge(2, 3);
    lifted_graph.insertEdge(4, 5);
    lifted_graph.insertEdge(5, 6);
    lifted_graph.insertEdge(6, 7);
    lifted_graph.insertEdge(0, 4);
    lifted_graph.insertEdge(1, 5);
    lifted_graph.insertEdge(2, 6);
    lifted_graph.insertEdge(3, 7);
    lifted_graph.insertEdge(0, 7);

    std::vector<double> edge_values(lifted_graph.numberOfEdges());
    edge_values[0] = -5;
    edge_values[1] = -5;
    edge_values[2] = -5;
    edge_values[3] = 5;
    edge_values[4] = .5;
    edge_values[5] = 5;
    edge_values[6] = 1;
    edge_values[7] = -1;
    edge_values[8] = 1;
    edge_values[9] = .5;
    edge_values[10] = -4;

    std::vector<char> edge_labels(lifted_graph.numberOfEdges());
    andres::graph::multicut_lifted::ilp_callback<andres::ilp::GurobiCallback>(original_graph, lifted_graph, edge_values, edge_labels, edge_labels);

    test(edge_labels[0] == 1);
    test(edge_labels[1] == 1);
    test(edge_labels[2] == 1);
    test(edge_labels[3] == 0);
    test(edge_labels[4] == 1);
    test(edge_labels[5] == 0);
    test(edge_labels[6] == 0);
    test(edge_labels[7] == 1);
    test(edge_labels[8] == 0);
    test(edge_labels[9] == 1);
    test(edge_labels[10] == 1);
}

int main() {
    test_3();
    test_4();
    test_8();
    return 0;
}
