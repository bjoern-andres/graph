#include <stdexcept>
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/cut-vertices.hxx"


inline void test(bool pred)
{ 
    if(!pred)
        throw std::runtime_error("Test failed."); 
}

using namespace andres::graph;

void test()
{
    Graph<> graph(8);
    graph.insertEdge(0, 1);
    graph.insertEdge(0, 2);
    graph.insertEdge(1, 3);
    graph.insertEdge(2, 4);
    graph.insertEdge(2, 5);
    graph.insertEdge(4, 6);
    graph.insertEdge(5, 6);

    CutVertices<decltype(graph)> cv(graph);
    cv.run();

    test(cv.isCutVertex(0) == true);
    test(cv.isCutVertex(1) == true);
    test(cv.isCutVertex(2) == true);
    test(cv.isCutVertex(3) == false);
    test(cv.isCutVertex(4) == false);
    test(cv.isCutVertex(5) == false);
    test(cv.isCutVertex(6) == false);
    test(cv.isCutVertex(7) == false);


    struct mask
    {
        bool vertex(std::size_t i) const
        {
            return !(i == 4 || i == 5 || i == 6);
        }

        bool edge(std::size_t i) const
        {
            return true;
        }
    };

    cv.run(mask());

    test(cv.isCutVertex(0) == true);
    test(cv.isCutVertex(1) == true);
    test(cv.isCutVertex(2) == false);
    test(cv.isCutVertex(3) == false);
    test(cv.isCutVertex(4) == false);
    test(cv.isCutVertex(5) == false);
    test(cv.isCutVertex(6) == false);
}

void testCompleteGraph()
{
    CompleteGraph<> graph(5);

    CutVertices<decltype(graph)> cv(graph);
    cv.run();

    test(cv.isCutVertex(0) == false);
    test(cv.isCutVertex(1) == false);
    test(cv.isCutVertex(2) == false);
    test(cv.isCutVertex(3) == false);
    test(cv.isCutVertex(4) == false);

    struct mask
    {
        bool vertex(std::size_t i) const
        {
            return true;
        }

        bool edge(std::size_t i) const
        {
            return !(i == 1 || i == 2 || i == 3 || i == 5 || i == 6);
        }
    };

    std::cout << graph.findEdge(1, 3).second << std::endl;
    std::cout << graph.findEdge(1, 4).second << std::endl;

    cv.run(mask());

    test(cv.isCutVertex(0) == false);
    test(cv.isCutVertex(1) == true);
    test(cv.isCutVertex(2) == true);
    test(cv.isCutVertex(3) == false);
    test(cv.isCutVertex(4) == false);
}

int main()
{
    test();

    testCompleteGraph();

    return 0;
}