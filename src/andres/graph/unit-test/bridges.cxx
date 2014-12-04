#include <stdexcept>
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/bridges.hxx"


inline void test(bool pred)
{ 
    if(!pred)
        throw std::runtime_error("Test failed."); 
}

using namespace andres::graph;

void test()
{
    Graph<> graph(8);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(0, 2); // 1
    graph.insertEdge(1, 3); // 2
    graph.insertEdge(2, 4); // 3
    graph.insertEdge(2, 5); // 4
    graph.insertEdge(4, 6); // 5
    graph.insertEdge(5, 6); // 6

    Bridges<decltype(graph)> br(graph);
    br.run();

    test(br.isBridge(0) == true);
    test(br.isBridge(1) == true);
    test(br.isBridge(2) == true);
    test(br.isBridge(3) == false);
    test(br.isBridge(4) == false);
    test(br.isBridge(5) == false);
    test(br.isBridge(6) == false);


    struct mask
    {
        bool vertex(std::size_t i) const
        {
            return true;
        }

        bool edge(std::size_t i) const
        {
            return !(i == 3);
        }
    };
    
    br.run(mask());

    test(br.isBridge(0) == true);
    test(br.isBridge(1) == true);
    test(br.isBridge(2) == true);
    test(br.isBridge(3) == false);
    test(br.isBridge(4) == true);
    test(br.isBridge(5) == true);
    test(br.isBridge(6) == true);
}

void testCompleteGraph()
{
    CompleteGraph<> graph(5);

    Bridges<decltype(graph)> br(graph);
    br.run();

    test(br.isBridge(graph.findEdge(0, 1).second) == false);
    test(br.isBridge(graph.findEdge(0, 2).second) == false);
    test(br.isBridge(graph.findEdge(0, 3).second) == false);
    test(br.isBridge(graph.findEdge(0, 4).second) == false);
    test(br.isBridge(graph.findEdge(1, 2).second) == false);
    test(br.isBridge(graph.findEdge(1, 3).second) == false);
    test(br.isBridge(graph.findEdge(1, 4).second) == false);
    test(br.isBridge(graph.findEdge(2, 3).second) == false);
    test(br.isBridge(graph.findEdge(2, 4).second) == false);
    test(br.isBridge(graph.findEdge(3, 4).second) == false);

    struct mask
    {
        bool vertex(std::size_t i) const
        {
            return true;
        }

        bool edge(std::size_t i) const
        {
            return !(i == 1 || i == 2 || i == 3);
        }
    };

    br.run(mask());

    test(br.isBridge(graph.findEdge(0, 1).second) == true);
    test(br.isBridge(graph.findEdge(0, 2).second) == false);
    test(br.isBridge(graph.findEdge(0, 3).second) == false);
    test(br.isBridge(graph.findEdge(0, 4).second) == false);
    test(br.isBridge(graph.findEdge(1, 2).second) == false);
    test(br.isBridge(graph.findEdge(1, 3).second) == false);
    test(br.isBridge(graph.findEdge(1, 4).second) == false);
    test(br.isBridge(graph.findEdge(2, 3).second) == false);
    test(br.isBridge(graph.findEdge(2, 4).second) == false);
    test(br.isBridge(graph.findEdge(3, 4).second) == false);
}

int main()
{
    test();

    testCompleteGraph();

    return 0;
}