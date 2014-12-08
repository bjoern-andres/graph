#include <stdexcept>
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/minimum-spanning-tree.hxx"


inline void test(bool pred)
{ 
    if(!pred)
        throw std::runtime_error("Test failed."); 
}

using namespace andres::graph;

void testMulticut()
{
    // a simple weighted graph in which an optimal multicut is non-trivial

    Graph<> graph;
    graph.insertVertices(6);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(0, 3); // 1
    graph.insertEdge(1, 2); // 2
    graph.insertEdge(1, 4); // 3
    graph.insertEdge(2, 5); // 4
    graph.insertEdge(3, 4); // 5
    graph.insertEdge(4, 5); // 6

    std::vector<int> weights(7);
    weights[0] = 5;
    weights[1] = -20;
    weights[2] = 5;
    weights[3] = 5;
    weights[4] = -20;
    weights[5] = 5;
    weights[6] = 5;

    std::vector<std::size_t> pred(graph.numberOfVertices());
    auto mst_value = findMSTSparseGraph(graph, weights, pred);

    test(mst_value == -25);

    struct mask
    {
        constexpr bool vertex(std::size_t i) const
        {
            return true;
        }

        bool edge(std::size_t i) const
        {
            return !(i == 1);
        }
    };

    mst_value = findMSTSparseGraph(graph, weights, mask(), pred);

    test(mst_value == 0);
}

void testMulticutCompleteGraph()
{
    CompleteGraph<> graph(5);
    
    std::vector<int> weights(10);
    weights[graph.findEdge(0, 1).second] = 10;
    weights[graph.findEdge(0, 2).second] = -1;
    weights[graph.findEdge(0, 3).second] = -1;
    weights[graph.findEdge(0, 4).second] = -1;
    weights[graph.findEdge(1, 2).second] = 10;
    weights[graph.findEdge(1, 3).second] = -1;
    weights[graph.findEdge(1, 4).second] = 4;
    weights[graph.findEdge(2, 3).second] = -1;
    weights[graph.findEdge(2, 4).second] = -1;
    weights[graph.findEdge(3, 4).second] = 10;

    std::vector<std::size_t> pred(graph.numberOfVertices());
    auto mst_value = findMSTDenseGraph(graph, weights, pred);

    test(mst_value == -4);

    struct mask
    {
        constexpr bool vertex(std::size_t i) const
        {
            return true;
        }

        bool edge(std::size_t i) const
        {
            return !(i == 5);
        }
    };

    mst_value = findMSTDenseGraph(graph, weights, mask(), pred);

    test(mst_value == 1);
}

int main()
{
    testMulticut();

    return 0;
}