#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/paths.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

struct SubgraphMask {
    bool vertex(const size_t v) const { return true; }
    bool edge(const size_t e) const { return e != 7; }
};

int main() {
    typedef andres::graph::Graph<> Graph;

    Graph graph(7);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(1, 2); // 1
    graph.insertEdge(2, 3); // 2
    graph.insertEdge(3, 4); // 3
    graph.insertEdge(4, 5); // 4
    graph.insertEdge(1, 6); // 5
    graph.insertEdge(6, 2); // 6
    graph.insertEdge(3, 5); // 7

    {
        size_t path[] =  {0, 1, 2, 3, 4, 5};

        std::pair<bool, size_t> p = andres::graph::findChord(graph, path, path + 6);
        test(p.first == true);
        test(p.second == 7);

        p = andres::graph::findChord(graph, path, path + 5);
        test(p.first == false);

        p = andres::graph::findChord(graph, path + 3, path + 6);
        test(p.first == true);
        test(p.second == 7);

        p = andres::graph::findChord(graph, SubgraphMask(), path + 3, path + 6, false);
        test(p.first == false);
    }
    {
        size_t path[] = {1, 6, 2};

        std::pair<bool, size_t> p = andres::graph::findChord(graph, path, path + 3);
        test(p.first == true);
        test(p.second == 1);

        p = andres::graph::findChord(graph, path, path + 3, true);
        test(p.first == false);
    }

    return 0;
}
