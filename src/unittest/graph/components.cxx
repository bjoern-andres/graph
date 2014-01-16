#include <cstddef>
#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/components.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

struct SubgraphMask {
    bool vertex(const std::size_t v) const { return v != 3; }
    bool edge(const std::size_t e) const { return e != 0; }
};

template<class COMPONENTS>
void testComponents() {
    typedef COMPONENTS Components;
    typedef typename Components::Graph Graph;

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
