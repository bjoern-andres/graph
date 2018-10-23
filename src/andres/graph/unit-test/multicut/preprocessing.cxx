#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/preprocessing.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

void testPreprocessing() {

    andres::graph::Graph<> graph;
    graph.insertVertices(8);
    graph.insertEdge(0, 1); // 0
    graph.insertEdge(0, 2); // 1
    graph.insertEdge(1, 2); // 2
    graph.insertEdge(1, 4); // 3
    graph.insertEdge(1, 3); // 4
    graph.insertEdge(2, 4); // 5
    graph.insertEdge(3, 4); // 6
    graph.insertEdge(3, 5); // 7
    graph.insertEdge(4, 5); // 8
    graph.insertEdge(1, 5); // 9
    graph.insertEdge(2, 6); // 10
    graph.insertEdge(2, 7); // 11
    graph.insertEdge(4, 6); // 12
    graph.insertEdge(4, 7); // 13
    graph.insertEdge(6, 7); // 14

    std::vector<double> weights(15);
    weights[0] = -2;
    weights[1] = 1;
    weights[2] = 2;
    weights[3] = 3;
    weights[4] = 3;
    weights[5] = -2;
    weights[6] = -2;
    weights[7] = 2;
    weights[8] = 3;
    weights[9] = 4;
    weights[10] = -2;
    weights[11] = -1;
    weights[12] = -3;
    weights[13] = -4;
    weights[14] = 2;

    std::cout << "Input graph:" << std::endl;
    std::cout << "|V| = " << graph.numberOfVertices() << ", Edges: " << std::endl;
    for (size_t e = 0; e < graph.numberOfEdges(); e++)
        std::cout << graph.vertexOfEdge(e, 0) << " - " << graph.vertexOfEdge(e, 1) << " : " << weights[e] << std::endl;

    std::cout << "-- Preprocessing -- " << std::endl;

    auto reduced_instance = andres::graph::multicut::preprocessing(graph, weights);

    auto constr = std::get<3>(reduced_instance);
    auto red_graph = std::get<0>(reduced_instance);
    auto red_costs = std::get<1>(reduced_instance);

    std::cout << "Contraints:" << std::endl;
    for (size_t i = 0; i < constr.size(); i++)
        std::cout << constr[i].first.first << " - " << constr[i].first.second << " : " << static_cast<int>(constr[i].second) << std::endl;

    std::cout << "Reduced graph:" << std::endl;
    std::cout << "|V| = " << red_graph.numberOfVertices() << ", Edges: " << std::endl;
    for (size_t e = 0; e < red_graph.numberOfEdges(); e++)
        std::cout << red_graph.vertexOfEdge(e, 0) << " - " << red_graph.vertexOfEdge(e, 1) << " : " << red_costs[e] << std::endl;    

}

int main()
{
    testPreprocessing();

    return 0;
}
