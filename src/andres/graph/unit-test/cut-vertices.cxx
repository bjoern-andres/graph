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

    std::vector<char> isCutVertex(graph.numberOfVertices());
    findCutVertices(graph, isCutVertex);

    test(isCutVertex[0] == true);
    test(isCutVertex[1] == true);
    test(isCutVertex[2] == true);
    test(isCutVertex[3] == false);
    test(isCutVertex[4] == false);
    test(isCutVertex[5] == false);
    test(isCutVertex[6] == false);
    test(isCutVertex[7] == false);


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

    findCutVertices(graph, mask(), isCutVertex);

    test(isCutVertex[0] == true);
    test(isCutVertex[1] == true);
    test(isCutVertex[2] == false);
    test(isCutVertex[3] == false);
    test(isCutVertex[4] == false);
    test(isCutVertex[5] == false);
    test(isCutVertex[6] == false);
}

void testCompleteGraph()
{
    CompleteGraph<> graph(5);

    std::vector<char> isCutVertex(graph.numberOfVertices());
    findCutVertices(graph, isCutVertex);

    test(isCutVertex[0] == false);
    test(isCutVertex[1] == false);
    test(isCutVertex[2] == false);
    test(isCutVertex[3] == false);
    test(isCutVertex[4] == false);

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

    findCutVertices(graph, mask(), isCutVertex);

    test(isCutVertex[0] == false);
    test(isCutVertex[1] == true);
    test(isCutVertex[2] == true);
    test(isCutVertex[3] == false);
    test(isCutVertex[4] == false);
}

int main()
{
    test();

    testCompleteGraph();

    return 0;
}