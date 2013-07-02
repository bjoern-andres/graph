#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/shortest-paths.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

struct SubgraphMask1 {
    bool vertex(const size_t v) const { return v != 1; }
    bool edge(const size_t e) const { return e != 3; }
};

struct SubgraphMask2 {
    bool vertex(const size_t v) const { return v != 1; }
    bool edge(const size_t e) const { return e < 3 || e > 6; }
};

struct SubgraphMask3 {
    bool vertex(const size_t v) const { return v != 1; }
    bool edge(const size_t e) const { return e != 3; }
};

struct SubgraphMask4 {
    bool vertex(const size_t j) const { return true; }
    bool edge(const size_t j) const { return j != 3; }
};

struct SubgraphMask5 {
    bool vertex(const size_t j) const { return true; }
    bool edge(const size_t j) const { return j != 3; }
};

struct SubgraphMask6 {
    bool vertex(const size_t j) const { return true; }
    bool edge(const size_t j) const { return j != 3; }
};

struct SubgraphMask7 {
    bool vertex(const size_t j) const { return true; }
    bool edge(const size_t j) const { return j != 3; }
};

int main() {
    // spsp, undirected, unweighted graph (breadth-first search)
    {
        andres::graph::Graph<> g(9);
        g.insertEdge(0, 1);
        g.insertEdge(1, 2);
        g.insertEdge(2, 3);
        g.insertEdge(0, 4);
        g.insertEdge(3, 4);
        g.insertEdge(0, 5);
        g.insertEdge(5, 3);
        g.insertEdge(0, 6);
        g.insertEdge(6, 7);
        g.insertEdge(7, 3);

        std::deque<size_t> path;

        bool found = andres::graph::spsp(g, 0, 3, path);
        test(found == true);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 4);
        test(path[2] == 3);

        found = andres::graph::spsp(g, 0, 8, path);
        test(found == false);
        test(path.size() == 0);

        found = andres::graph::spsp(g, 0, 0, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 0);

        found = andres::graph::spsp(g, 8, 8, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 8);
        
        found = andres::graph::spspEdges(g, 0, 3, path);
        test(found == true);
        test(path.size() == 2);
        test(path[0] == 3);
        test(path[1] == 4);
        
        found = andres::graph::spspEdges(g, 0, 8, path);
        test(found == false);
        test(path.size() == 0);
        
        found = andres::graph::spspEdges(g, 0, 0, path);
        test(found == true);
        test(path.size() == 0);
        
        found = andres::graph::spspEdges(g, 8, 8, path);
        test(found == true);
        test(path.size() == 0);
    }

    // spsp, undirected, weighted graph (Dijkstra)
    {
        andres::graph::Graph<> g(9);
        g.insertEdge(0, 1); // 0
        g.insertEdge(1, 2); // 1
        g.insertEdge(2, 3); // 2
        g.insertEdge(0, 4); // 3
        g.insertEdge(3, 4); // 4
        g.insertEdge(0, 5); // 5
        g.insertEdge(5, 3); // 6
        g.insertEdge(0, 6); // 7
        g.insertEdge(6, 7); // 8
        g.insertEdge(7, 3); // 9

        std::vector<float> edgeWeights(10);
        edgeWeights[0] = 0.1f;
        edgeWeights[1] = 0.2f;
        edgeWeights[2] = 0.3f;
        edgeWeights[3] = 0.1f;
        edgeWeights[4] = 0.6f;
        edgeWeights[5] = 0.3f;
        edgeWeights[6] = 0.4f;
        edgeWeights[7] = 0.1f;
        edgeWeights[8] = 0.2f;
        edgeWeights[9] = 1.0f;

        std::deque<size_t> path;
        float distance = 0;      

        // start vertex 0
        andres::graph::spsp(g, 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 0);

        andres::graph::spsp(g, 0, 1, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 1);
        
        andres::graph::spsp(g, 0, 2, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 1);
        test(path[2] == 2);
        
        andres::graph::spsp(g, 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f + 0.3f);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 1);
        test(path[2] == 2);
        test(path[3] == 3);
        
        andres::graph::spsp(g, 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 4);
        
        andres::graph::spsp(g, 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 5);
        
        andres::graph::spsp(g, 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 6);
        
        andres::graph::spsp(g, 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 6);
        test(path[2] == 7);
        
        andres::graph::spsp(g, 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 0, 1, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 0);
        
        andres::graph::spspEdges(g, 0, 2, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 1);
        
        andres::graph::spspEdges(g, 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f + 0.3f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 1);
        test(path[2] == 2);
        
        andres::graph::spspEdges(g, 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 3);
        
        andres::graph::spspEdges(g, 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 5);
        
        andres::graph::spspEdges(g, 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 7);
        
        andres::graph::spspEdges(g, 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 2);
        test(path[0] == 7);
        test(path[1] == 8);
        
        andres::graph::spspEdges(g, 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        // start vertex 2
        andres::graph::spsp(g, 2, 0, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f);
        test(path.size() == 3);
        test(path[0] == 2);
        test(path[1] == 1);
        test(path[2] == 0);

        andres::graph::spsp(g, 2, 1, edgeWeights.begin(), path, distance);
        test(distance == 0.2f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 1);
        
        andres::graph::spsp(g, 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spsp(g, 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 3);
        
        andres::graph::spsp(g, 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.1f);
        test(path.size() == 4);
        test(path[0] == 2);
        test(path[1] == 1);
        test(path[2] == 0);
        test(path[3] == 4);
        
        andres::graph::spsp(g, 2, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.3f);
        test(path.size() == 4);
        test(path[0] == 2);
        test(path[1] == 1);
        test(path[2] == 0);
        test(path[3] == 5);
        
        andres::graph::spsp(g, 2, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.1f);
        test(path.size() == 4);
        test(path[0] == 2);
        test(path[1] == 1);
        test(path[2] == 0);
        test(path[3] == 6);
        
        andres::graph::spsp(g, 2, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.1f + 0.2f);
        test(path.size() == 5);
        test(path[0] == 2);
        test(path[1] == 1);
        test(path[2] == 0);
        test(path[3] == 6);
        test(path[4] == 7);
        
        andres::graph::spsp(g, 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 0, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f);
        test(path.size() == 2);
        test(path[0] == 1);
        test(path[1] == 0);
        
        andres::graph::spspEdges(g, 2, 1, edgeWeights.begin(), path, distance);
        test(distance == 0.2f);
        test(path.size() == 1);
        test(path[0] == 1);
        
        andres::graph::spspEdges(g, 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spspEdges(g, 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.1f);
        test(path.size() == 3);
        test(path[0] == 1);
        test(path[1] == 0);
        test(path[2] == 3);
        
        andres::graph::spspEdges(g, 2, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.3f);
        test(path.size() == 3);
        test(path[0] == 1);
        test(path[1] == 0);
        test(path[2] == 5);
        
        andres::graph::spspEdges(g, 2, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.1f);
        test(path.size() == 3);
        test(path[0] == 1);
        test(path[1] == 0);
        test(path[2] == 7);
        
        andres::graph::spspEdges(g, 2, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.2f + 0.1f + 0.1f + 0.2f);
        test(path.size() == 4);
        test(path[0] == 1);
        test(path[1] == 0);
        test(path[2] == 7);
        test(path[3] == 8);
        
        andres::graph::spspEdges(g, 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        // with subgraph mask
        // start vertex 0
        andres::graph::spsp(g, SubgraphMask1(), 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 0);

        andres::graph::spsp(g, SubgraphMask1(), 0, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 2, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.3f);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);
        test(path[3] == 2);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.6f);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);
        test(path[3] == 4);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 5);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 6);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 6);
        test(path[2] == 7);
        
        andres::graph::spsp(g, SubgraphMask1(), 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 2, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.3f);
        test(path.size() == 3);
        test(path[0] == 5);
        test(path[1] == 6);
        test(path[2] == 2);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f);
        test(path.size() == 2);
        test(path[0] == 5);
        test(path[1] == 6);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.6f);
        test(path.size() == 3);
        test(path[0] == 5);
        test(path[1] == 6);
        test(path[2] == 4);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 5);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 7);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 2);
        test(path[0] == 7);
        test(path[1] == 8);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        // start vertex 2
        andres::graph::spsp(g, SubgraphMask1(), 2, 0, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.3f);
        test(path.size() == 4);
        test(path[0] == 2);
        test(path[1] == 3);
        test(path[2] == 5);
        test(path[3] == 0);

        andres::graph::spsp(g, SubgraphMask1(), 2, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask1(), 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spsp(g, SubgraphMask1(), 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 3);
        
        andres::graph::spsp(g, SubgraphMask1(), 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.6f);
        test(path.size() == 3);
        test(path[0] == 2);
        test(path[1] == 3);
        test(path[2] == 4);
        
        andres::graph::spsp(g, SubgraphMask1(), 2, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f);
        test(path.size() == 3);
        test(path[0] == 2);
        test(path[1] == 3);
        test(path[2] == 5);
        
        andres::graph::spsp(g, SubgraphMask1(), 2, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.3f + 0.1f);
        test(path.size() == 5);
        test(path[0] == 2);
        test(path[1] == 3);
        test(path[2] == 5);
        test(path[3] == 0);
        test(path[4] == 6);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 0, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.3f);
        test(path.size() == 3);
        test(path[0] == 2);
        test(path[1] == 6);
        test(path[2] == 5);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.6f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 4);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 6);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.3f + 0.1f);
        test(path.size() == 4);
        test(path[0] == 2);
        test(path[1] == 6);
        test(path[2] == 5);
        test(path[3] == 7);
        
        // 2 to 7 is ambiguous. therefore not part of the unit test

        andres::graph::spsp(g, SubgraphMask1(), 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask1(), 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
    }

    // spsp, directed, unweighted graph (breadth-first search)
    {
        andres::graph::Digraph<> g(9);
        g.insertEdge(0, 1); // 0
        g.insertEdge(1, 2); // 1
        g.insertEdge(2, 3); // 2
        g.insertEdge(0, 4); // 3
        g.insertEdge(3, 4); // 4
        g.insertEdge(0, 5); // 5
        g.insertEdge(5, 3); // 6
        g.insertEdge(0, 6); // 7
        g.insertEdge(6, 7); // 8
        g.insertEdge(7, 3); // 9

        std::deque<size_t> path;

        bool found = andres::graph::spsp(g, 0, 3, path);   
        test(found == true);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);

        found = andres::graph::spsp(g, 3, 0, path);
        test(found == false);
        test(path.size() == 0);

        found = andres::graph::spsp(g, 0, 8, path);
        test(found == false);
        test(path.size() == 0);

        found = andres::graph::spsp(g, 0, 0, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 0);

        found = andres::graph::spsp(g, 8, 8, path);
        test(found == true);
        test(path.size() == 1);
        test(path[0] == 8);
        
        found = andres::graph::spspEdges(g, 0, 3, path);
        test(found == true);
        test(path.size() == 2);
        test(path[0] == 5);
        test(path[1] == 6);
        
        found = andres::graph::spspEdges(g, 3, 0, path);
        test(found == false);
        test(path.size() == 0);
        
        found = andres::graph::spspEdges(g, 0, 8, path);
        test(found == false);
        test(path.size() == 0);
        
        found = andres::graph::spspEdges(g, 0, 0, path);
        test(found == true);
        test(path.size() == 0);
        
        found = andres::graph::spspEdges(g, 8, 8, path);
        test(found == true);
        test(path.size() == 0);

        // with subgraph mask

        found = andres::graph::spsp(g, SubgraphMask2(), 0, 3, path);
        test(found == true);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 6);
        test(path[2] == 7);
        test(path[3] == 3);
        
        found = andres::graph::spspEdges(g, SubgraphMask2(), 0, 3, path);
        test(found == true);
        test(path.size() == 3);
        test(path[0] == 7);
        test(path[1] == 8);
        test(path[2] == 9);
    }

    // spsp, directed, weighted graph (Dijkstra)
    {
        andres::graph::Digraph<> g(9);
        g.insertEdge(0, 1); // 0
        g.insertEdge(1, 2); // 1
        g.insertEdge(2, 3); // 2
        g.insertEdge(0, 4); // 3
        g.insertEdge(3, 4); // 4
        g.insertEdge(0, 5); // 5
        g.insertEdge(5, 3); // 6
        g.insertEdge(0, 6); // 7
        g.insertEdge(6, 7); // 8
        g.insertEdge(7, 3); // 9

        std::vector<float> edgeWeights(10);
        edgeWeights[0] = 0.1f;
        edgeWeights[1] = 0.2f;
        edgeWeights[2] = 0.3f;
        edgeWeights[3] = 0.1f;
        edgeWeights[4] = 0.6f;
        edgeWeights[5] = 0.3f;
        edgeWeights[6] = 0.4f;
        edgeWeights[7] = 0.1f;
        edgeWeights[8] = 0.2f;
        edgeWeights[9] = 1.0f;

        std::deque<size_t> path;
        float distance = 0;      

        // start vertex 0
        andres::graph::spsp(g, 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 0);

        andres::graph::spsp(g, 0, 1, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 1);
        
        andres::graph::spsp(g, 0, 2, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 1);
        test(path[2] == 2);
        
        andres::graph::spsp(g, 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f + 0.3f);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 1);
        test(path[2] == 2);
        test(path[3] == 3);
        
        andres::graph::spsp(g, 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 4);
        
        andres::graph::spsp(g, 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 5);
        
        andres::graph::spsp(g, 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 6);
        
        andres::graph::spsp(g, 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 6);
        test(path[2] == 7);
        
        andres::graph::spsp(g, 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 0, 1, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 0);
        
        andres::graph::spspEdges(g, 0, 2, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 1);
        
        andres::graph::spspEdges(g, 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f + 0.3f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 1);
        test(path[2] == 2);
        
        andres::graph::spspEdges(g, 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 3);
        
        andres::graph::spspEdges(g, 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 5);
        
        andres::graph::spspEdges(g, 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 7);
        
        andres::graph::spspEdges(g, 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 2);
        test(path[0] == 7);
        test(path[1] == 8);
        
        andres::graph::spspEdges(g, 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        // start vertex 2
        andres::graph::spsp(g, 2, 0, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        andres::graph::spsp(g, 2, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spsp(g, 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 3);
        
        andres::graph::spsp(g, 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.6f);
        test(path.size() == 3);
        test(path[0] == 2);
        test(path[1] == 3);
        test(path[2] == 4);
        
        andres::graph::spsp(g, 2, 5, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, 2, 6, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, 2, 7, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        andres::graph::spspEdges(g, 2, 0, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spspEdges(g, 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.6f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 4);
        
        andres::graph::spspEdges(g, 2, 5, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 6, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 7, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        // with subgraph mask
        // start vertex 0
        andres::graph::spsp(g, SubgraphMask3(), 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 0);

        andres::graph::spsp(g, SubgraphMask3(), 0, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 2, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.6f);
        test(path.size() == 4);
        test(path[0] == 0);
        test(path[1] == 5);
        test(path[2] == 3);
        test(path[3] == 4);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 5);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 2);
        test(path[0] == 0);
        test(path[1] == 6);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 3);
        test(path[0] == 0);
        test(path[1] == 6);
        test(path[2] == 7);
        
        andres::graph::spsp(g, SubgraphMask3(), 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 0, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 2, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f);
        test(path.size() == 2);
        test(path[0] == 5);
        test(path[1] == 6);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.4f + 0.6f);
        test(path.size() == 3);
        test(path[0] == 5);
        test(path[1] == 6);
        test(path[2] == 4);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 5, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 5);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 6, edgeWeights.begin(), path, distance);
        test(distance == 0.1f);
        test(path.size() == 1);
        test(path[0] == 7);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 7, edgeWeights.begin(), path, distance);
        test(distance == 0.1f + 0.2f);
        test(path.size() == 2);
        test(path[0] == 7);
        test(path[1] == 8);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 0, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        // start vertex 2
        andres::graph::spsp(g, SubgraphMask3(), 2, 0, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);

        andres::graph::spsp(g, SubgraphMask3(), 2, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 3);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.6f);
        test(path.size() == 3);
        test(path[0] == 2);
        test(path[1] == 3);
        test(path[2] == 4);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 5, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 6, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 7, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spsp(g, SubgraphMask3(), 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 0, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 1, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 2, edgeWeights.begin(), path, distance);
        test(distance == 0);
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 3, edgeWeights.begin(), path, distance);
        test(distance == 0.3f);
        test(path.size() == 1);
        test(path[0] == 2);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 4, edgeWeights.begin(), path, distance);
        test(distance == 0.3f + 0.6f);
        test(path.size() == 2);
        test(path[0] == 2);
        test(path[1] == 4);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 5, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 6, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 7, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
        
        andres::graph::spspEdges(g, SubgraphMask3(), 2, 8, edgeWeights.begin(), path, distance);
        test(distance == std::numeric_limits<float>::infinity());
        test(path.size() == 0);
    }

    // sssp, undirected graph (Dijkstra)
    {
        andres::graph::Graph<> g(5);
        g.insertEdge(0, 1); // 0
        g.insertEdge(1, 2); // 1
        g.insertEdge(2, 3); // 2
        g.insertEdge(1, 4); // 3
        g.insertEdge(4, 3); // 4

        // unweighted graph
        {
            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, 0, distances.begin(), parents.begin());

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 2);

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 1);

            // convenience function without parents
            andres::graph::sssp(g, 0, distances.begin());

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 2);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, 0, distances.begin(), parents.begin(), parentsEdges.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 2);
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 1);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 2);
            test(parentsEdges[4] == 3);
            
            // convenience function without parents
            andres::graph::ssspEdges(g, 0, distances.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 2);
        }

        // unweighted subgraph
        {
            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, SubgraphMask4(), 0, distances.begin(), parents.begin());

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 4);

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 3);

            // convenience function without parents
            andres::graph::sssp(g, SubgraphMask4(), 0, distances.begin());

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 4);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, SubgraphMask4(), 0, distances.begin(), parents.begin(), parentsEdges.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 4);
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 3);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 2);
            test(parentsEdges[4] == 4);
            
            
            // convenience function without parents
            andres::graph::ssspEdges(g, SubgraphMask4(), 0, distances.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 4);
        }

        std::vector<unsigned int> edgeWeights(5);
        edgeWeights[0] = 1;
        edgeWeights[1] = 1;
        edgeWeights[2] = 3;
        edgeWeights[3] = 2;
        edgeWeights[4] = 1;

        // weighted graph
        {
            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, 0, edgeWeights.begin(), distances.begin(), 
                parents.begin()
            );

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 4);
            test(distances[4] == 3);

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 4);
            test(parents[4] == 1);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, 0, edgeWeights.begin(), distances.begin(), parents.begin(), parentsEdges.begin()
            );
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 4);
            test(distances[4] == 3);
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 4);
            test(parents[4] == 1);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 4);
            test(parentsEdges[4] == 3);
        }

        // weighted subgraph
        {


            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, SubgraphMask5(), 0, edgeWeights.begin(),
                distances.begin(), parents.begin()
            );

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 5);
            test(distances[4] == 6);

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 3);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, SubgraphMask5(), 0, edgeWeights.begin(), distances.begin(), parents.begin(), parentsEdges.begin()
            );
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 5);
            test(distances[4] == 6);
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 3);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 2);
            test(parentsEdges[4] == 4);
        }
    }

    // sssp, directed graph (Dijkstra)
    {
        andres::graph::Digraph<> g(5);
        g.insertEdge(0, 1); // 0
        g.insertEdge(1, 2); // 1
        g.insertEdge(2, 3); // 2
        g.insertEdge(1, 4); // 3
        g.insertEdge(4, 3); // 4
        g.insertEdge(4, 0); // 5

        // unweighted graph
        {
            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, 0, distances.begin(), parents.begin());

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 2);

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 1);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, 0, distances.begin(), parents.begin(), parentsEdges.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == 2);
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 1);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 2);
            test(parentsEdges[4] == 3);
        }

        // unweighted subgraph
        {


            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, SubgraphMask6(), 0, distances.begin(), parents.begin());

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == std::numeric_limits<unsigned int>::max());

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 0);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, SubgraphMask6(), 0, distances.begin(), parents.begin(), parentsEdges.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 3);
            test(distances[4] == std::numeric_limits<unsigned int>::max());
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 0);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 2);
            test(parentsEdges[4] == 0);
            
        }

        std::vector<unsigned int> edgeWeights(6);
        edgeWeights[0] = 1;
        edgeWeights[1] = 1;
        edgeWeights[2] = 3;
        edgeWeights[3] = 2;
        edgeWeights[4] = 1;
        edgeWeights[5] = 1;

        // weighted graph
        {
            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, 0, edgeWeights.begin(), distances.begin(), 
                parents.begin()
            );

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 4);
            test(distances[4] == 3);

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 4);
            test(parents[4] == 1);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, 0, edgeWeights.begin(), distances.begin(),
                                     parents.begin(), parentsEdges.begin()
                                     );
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 4);
            test(distances[4] == 3);
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 4);
            test(parents[4] == 1);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 4);
            test(parentsEdges[4] == 3);
        }

        // weighted subgraph
        {


            std::vector<unsigned int> distances(g.numberOfVertices());
            std::vector<size_t> parents(g.numberOfVertices());
            andres::graph::sssp(g, SubgraphMask7(), 0, edgeWeights.begin(),
                distances.begin(), parents.begin()
            );

            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 5);
            test(distances[4] == std::numeric_limits<unsigned int>::max());

            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 0);
            
            // edge output
            
            std::vector<size_t> parentsEdges(g.numberOfVertices());
            andres::graph::ssspEdges(g, SubgraphMask7(), 0, edgeWeights.begin(),distances.begin(), parents.begin(), parentsEdges.begin());
            
            test(distances[0] == 0);
            test(distances[1] == 1);
            test(distances[2] == 2);
            test(distances[3] == 5);
            test(distances[4] == std::numeric_limits<unsigned int>::max());
            
            test(parents[0] == 0);
            test(parents[1] == 0);
            test(parents[2] == 1);
            test(parents[3] == 2);
            test(parents[4] == 0);
            
            test(parentsEdges[0] == 0);
            test(parentsEdges[1] == 0);
            test(parentsEdges[2] == 1);
            test(parentsEdges[3] == 2);
            test(parentsEdges[4] == 0);
        }
    }

    return 0;
}
