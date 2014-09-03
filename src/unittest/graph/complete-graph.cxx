#include <stdexcept>

#include "andres/graph/complete-graph.hxx"

inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

typedef andres::graph::CompleteGraph<> CompleteGraph;
typedef andres::graph::Adjacency<> Adjacency;
typedef std::pair<bool, std::size_t> Pair;

void testNumberOfVerticesEdges() {
    for(std::size_t j = 0; j < 20; ++j) {
        CompleteGraph g(j);
        test(g.numberOfVertices() == j);
        test(g.numberOfEdges() == j * (j - 1) / 2);
        for(std::size_t k = 0; k < j; ++k) {
            test(g.numberOfEdgesFromVertex(k) == j - 1);
            test(g.numberOfEdgesToVertex(k) == j - 1);
        }
    }
}

void testFindEdge() {
    CompleteGraph g(4);
    { Pair p = g.findEdge(0, 0); test(p.first == false); }
    { Pair p = g.findEdge(0, 1); test(p.first == true); test(p.second == 0); }
    { Pair p = g.findEdge(0, 2); test(p.first == true); test(p.second == 1); }
    { Pair p = g.findEdge(0, 3); test(p.first == true); test(p.second == 2); }
    { Pair p = g.findEdge(1, 0); test(p.first == true); test(p.second == 0); }
    { Pair p = g.findEdge(1, 1); test(p.first == false); }
    { Pair p = g.findEdge(1, 2); test(p.first == true); test(p.second == 3); }
    { Pair p = g.findEdge(1, 3); test(p.first == true); test(p.second == 4); }
    { Pair p = g.findEdge(2, 0); test(p.first == true); test(p.second == 1); }
    { Pair p = g.findEdge(2, 1); test(p.first == true); test(p.second == 3); }
    { Pair p = g.findEdge(2, 2); test(p.first == false); }
    { Pair p = g.findEdge(2, 3); test(p.first == true); test(p.second == 5); }
    { Pair p = g.findEdge(3, 0); test(p.first == true); test(p.second == 2); }
    { Pair p = g.findEdge(3, 1); test(p.first == true); test(p.second == 4); }
    { Pair p = g.findEdge(3, 2); test(p.first == true); test(p.second == 5); }
    { Pair p = g.findEdge(3, 3); test(p.first == false); }
}

void testAdjacency() {
    CompleteGraph g(4);
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
            const Adjacency adjacency = g.adjacencyFromVertex(v, j);
            if(j < v) {
                test(adjacency.vertex() == j);
                test(adjacency.edge() == g.findEdge(v, adjacency.vertex()).second);
            }
            else {
                test(adjacency.vertex() == j + 1);
                test(adjacency.edge() == g.findEdge(v, adjacency.vertex()).second);
            }

            test(g.vertexFromVertex(v, j) == adjacency.vertex());
            test(g.edgeFromVertex(v, j) == adjacency.edge());
        }
        for(std::size_t j = 0; j < g.numberOfEdgesToVertex(v); ++j) {
            const Adjacency adjacency = g.adjacencyToVertex(v, j);
            if(j < v) {
                test(adjacency.vertex() == j);
                test(adjacency.edge() == g.findEdge(v, adjacency.vertex()).second);
            }
            else {
                test(adjacency.vertex() == j + 1);
                test(adjacency.edge() == g.findEdge(v, adjacency.vertex()).second);
            }

            test(g.vertexToVertex(v, j) == adjacency.vertex());
            test(g.edgeToVertex(v, j) == adjacency.edge());
        }
    }
}

void testVertexOfEdge() {
    CompleteGraph g(4);
    for(std::size_t e = 0; e < g.numberOfEdges(); ++ e) {
        const std::size_t v0 = g.vertexOfEdge(e, 0);
        const std::size_t v1 = g.vertexOfEdge(e, 1);
        test(v0 < v1);
        test(g.findEdge(v0, v1).second == e);
    }
}

int main() {
    testNumberOfVerticesEdges();
    testFindEdge();
    testAdjacency();
    testVertexOfEdge();

    return 0;
}
