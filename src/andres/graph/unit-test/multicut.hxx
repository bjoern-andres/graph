#pragma once

#include "andres/graph/graph.hxx"
#include "andres/graph/graph-complete.hxx"
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

template<class ILP>
void testMulticutCompleteGraph() {
    typedef ILP Ilp;
    typedef andres::graph::CompleteGraph<> Graph;
    typedef andres::graph::Multicut<Graph, Ilp> Multicut;

    // a simple weighted complete graph in which an optimal multicut is non-trivial

    Graph graph(5);

    std::vector<double> weights(10);
    weights[graph.findEdge(0, 1).second] = 10;
    weights[graph.findEdge(0, 2).second] = 1;
    weights[graph.findEdge(0, 3).second] = 1;
    weights[graph.findEdge(0, 4).second] = 1;
    weights[graph.findEdge(1, 2).second] = 10;
    weights[graph.findEdge(1, 3).second] = 1;
    weights[graph.findEdge(1, 4).second] = 1;
    weights[graph.findEdge(2, 3).second] = 1;
    weights[graph.findEdge(2, 4).second] = 1;
    weights[graph.findEdge(3, 4).second] = 10;

    Multicut mc;
    mc.ilp().setVerbosity(false);
    mc.setup(graph, weights);
    mc.solve(10);

    test(mc.label(graph.findEdge(0, 1).second) == 0);
    test(mc.label(graph.findEdge(0, 2).second) == 0);
    test(mc.label(graph.findEdge(0, 3).second) == 1);
    test(mc.label(graph.findEdge(0, 4).second) == 1);
    test(mc.label(graph.findEdge(1, 2).second) == 0);
    test(mc.label(graph.findEdge(1, 3).second) == 1);
    test(mc.label(graph.findEdge(1, 4).second) == 1);
    test(mc.label(graph.findEdge(2, 3).second) == 1);
    test(mc.label(graph.findEdge(2, 4).second) == 1);
    test(mc.label(graph.findEdge(3, 4).second) == 0);
}
