#include <stdexcept>
#include "andres/graph/graph.hxx"
#include "andres/graph/triangles.hxx"


inline void test(bool pred)
{ 
    if(!pred)
        throw std::runtime_error("Test failed."); 
}

using namespace andres::graph;

void testTriangles()
{
    andres::graph::Graph<> test_graph;
    test_graph.insertVertices(6);
    test_graph.insertEdge(0, 1);
    test_graph.insertEdge(0, 2);
    test_graph.insertEdge(1, 2);
    test_graph.insertEdge(1, 4);
    test_graph.insertEdge(1, 3);
    test_graph.insertEdge(3, 4);
    test_graph.insertEdge(2, 4);
    test_graph.insertEdge(4, 5);
    test_graph.insertEdge(3, 5);

    std::cout << "-- Test graph --" << std::endl;
    std::cout << "Number of vertices: " << test_graph.numberOfVertices() << std::endl;
    std::cout << "Edges:" << std::endl;
    for (size_t e = 0; e < test_graph.numberOfEdges(); e++)
    {
        auto v0 = test_graph.vertexOfEdge(e, 0);
        auto v1 = test_graph.vertexOfEdge(e, 1);
        std::cout << v0 << " - " << v1 << std::endl;
    }
    
    auto triangles = findTriangles(test_graph);

    test(triangles.size() == 4);

    std::cout << "Triangles found: " << std::endl;

    for (auto it = triangles.begin(); it != triangles.end(); it++)
    {
        auto triangle = *it;
        std::cout << triangle[0] << " - " << triangle[1] << " - " << triangle[2] << std::endl;
    }
}

int main()
{
    testTriangles();

    return 0;
}