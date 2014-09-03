#include <stdexcept>

#include "andres/graph/complete-graph.hxx"

inline void test(const bool& pred) {
    if(!pred) throw std::runtime_error("Test failed.");
}

typedef andres::graph::CompleteGraph<> CompleteGraph;
typedef CompleteGraph::VertexIterator VertexIterator;
typedef CompleteGraph::EdgeIterator EdgeIterator;
typedef CompleteGraph::AdjacencyIterator AdjacencyIterator;
typedef andres::graph::Adjacency<> Adjacency;
typedef std::pair<bool, std::size_t> Pair;

void testConstructionAndNumbers() {
    {
        const CompleteGraph g;
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }
    for(std::size_t j = 0; j < 20; ++j) {
        const CompleteGraph g(j);
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == j);
        test(g.numberOfEdges() == j * (j - 1) / 2);
        for(std::size_t k = 0; k < j; ++k) {
            test(g.numberOfEdgesFromVertex(k) == j - 1);
            test(g.numberOfEdgesToVertex(k) == j - 1);
        }
    }
}

// tests findEdge explicitly for the complete graph K4
void testFindEdge() {
    const CompleteGraph g(4);
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

// tests consistency of adjacency functions with findEdge
void testAdjacency() {
    const CompleteGraph g(4);
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
    for(std::size_t e = 0; e < g.numberOfEdges(); ++ e) {
        const std::size_t v0 = g.vertexOfEdge(e, 0);
        const std::size_t v1 = g.vertexOfEdge(e, 1);
        test(v0 < v1);
        test(g.findEdge(v0, v1).second == e);
    }
}

void testVertexIterator() {
    const CompleteGraph g(5);

    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // verticesFromVertex
        {
            // operator[]
            VertexIterator it = g.verticesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(it[j] == g.vertexFromVertex(v, j));
            }
        }
        {
            // operator*
            VertexIterator it = g.verticesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(*it == g.vertexFromVertex(v, j));
                ++it;
            }
            test(it == g.verticesFromVertexEnd(v));
        }
        // verticesToVertex
        {
            // operator[]
            VertexIterator it = g.verticesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(it[j] == g.vertexToVertex(v, j));
            }
        }
        {
            // operator*, operator++, operator==
            VertexIterator it = g.verticesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(*it == g.vertexToVertex(v, j));
                ++it;
            }
            test(it == g.verticesToVertexEnd(v));
        }
    }

    // operator==, operator!=
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        VertexIterator it = g.verticesFromVertexBegin(v);
        for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
            for(std::size_t w = 0; w < g.numberOfVertices(); ++w) {
                VertexIterator it2 = g.verticesFromVertexBegin(w);
                for(std::size_t k = 0; k < g.numberOfVertices() - 1; ++k) {
                    test((it + j == it2 + k) == (v == w && j == k));
                    test((it + j != it2 + k) == (v != w || j != k));
                }
            }
        }
    }

    // increment and decrement operators (consistency with operator==, operator!=)
    {
        VertexIterator it = g.verticesFromVertexBegin(1);
        // operator++ (prefix)
        VertexIterator it2 = ++it;
        test(it2 == it);
        test(*it == 2);
        // operator-- (prefix)
        it2 = --it;
        test(it2 == it);
        test(*it == 0);
        // operator++ (postfix)
        it2 = it++;
        test(it2 != it);
        test(*it == 2);
        test(*it2 == 0);
        // operator-- (postfix)
        it2 = it--;
        test(it2 != it);
        test(*it == 0);
        test(*it2 == 2);
        // operator+=
        it2 = (it += 2);
        test(it2 == it);
        test(*it == 3);
        // operator-=
        it2 = (it -= 2);
        test(it2 == it);
        test(*it == 0);
    }

    // operator+, operator-
    {
        VertexIterator it = g.verticesFromVertexBegin(1);
        // operator+
        VertexIterator it2 = it + 2;
        test(it2 != it);
        test(*it == 0);
        test(*it2 == 3);
        // operator-
        VertexIterator it3 = it2 - 2;
        test(it3 != it2);
        test(*it2 == 3);
        test(*it3 == 0);
    }

    // operator>, operator>=, operator<, operator<=
    {
        VertexIterator it = g.verticesFromVertexBegin(1);
        VertexIterator it2 = it;
        ++it2;
        test(!(it2 < it));
        test(!(it2 <= it));
        test(it2 > it);
        test(it2 >= it);
    }
}

void testEdgeIterator() {
    const CompleteGraph g(5);

    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // edgesFromVertex
        {
            // operator[]
            EdgeIterator it = g.edgesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(it[j] == g.edgeFromVertex(v, j));
            }
        }
        {
            // operator*
            EdgeIterator it = g.edgesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(*it == g.edgeFromVertex(v, j));
                ++it;
            }
            test(it == g.edgesFromVertexEnd(v));
        }
        // edgesToVertex
        {
            // operator[]
            EdgeIterator it = g.edgesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(it[j] == g.edgeToVertex(v, j));
            }
        }
        {
            // operator*, operator++, operator==
            EdgeIterator it = g.edgesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(*it == g.edgeToVertex(v, j));
                ++it;
            }
            test(it == g.edgesToVertexEnd(v));
        }
    }

    // operator==, operator!=
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        EdgeIterator it = g.edgesFromVertexBegin(v);
        for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
            for(std::size_t w = 0; w < g.numberOfVertices(); ++w) {
                EdgeIterator it2 = g.edgesFromVertexBegin(w);
                for(std::size_t k = 0; k < g.numberOfVertices() - 1; ++k) {
                    test((it + j == it2 + k) == (v == w && j == k));
                    test((it + j != it2 + k) == (v != w || j != k));
                }
            }
        }
    }

    // increment and decrement operators (consistency with operator==, operator!=)
    {
        EdgeIterator it = g.edgesFromVertexBegin(1);
        // operator++ (prefix)
        EdgeIterator it2 = ++it;
        test(it2 == it);
        test(*it == g.findEdge(1, 2).second);
        // operator-- (prefix)
        it2 = --it;
        test(it2 == it);
        test(*it == g.findEdge(1, 0).second);
        // operator++ (postfix)
        it2 = it++;
        test(it2 != it);
        test(*it == g.findEdge(1, 2).second);
        test(*it2 == g.findEdge(1, 0).second);
        // operator-- (postfix)
        it2 = it--;
        test(it2 != it);
        test(*it == g.findEdge(1, 0).second);
        test(*it2 == g.findEdge(1, 2).second);
        // operator+=
        it2 = (it += 2);
        test(it2 == it);
        test(*it == g.findEdge(1, 3).second);
        // operator-=
        it2 = (it -= 2);
        test(it2 == it);
        test(*it == g.findEdge(1, 0).second);
    }

    // operator+, operator-
    {
        EdgeIterator it = g.edgesFromVertexBegin(1);
        // operator+
        EdgeIterator it2 = it + 2;
        test(it2 != it);
        test(*it == g.findEdge(1, 0).second);
        test(*it2 == g.findEdge(1, 3).second);
        // operator-
        EdgeIterator it3 = it2 - 2;
        test(it3 != it2);
        test(*it2 == g.findEdge(1, 3).second);
        test(*it3 == g.findEdge(1, 0).second);
    }

    // operator>, operator>=, operator<, operator<=
    {
        EdgeIterator it = g.edgesFromVertexBegin(1);
        EdgeIterator it2 = it;
        ++it2;
        test(!(it2 < it));
        test(!(it2 <= it));
        test(it2 > it);
        test(it2 >= it);
    }
}

void testAdjacencyIterator() {
    const CompleteGraph g(5);

    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // adjacenciesFromVertex
        {
            // operator[]
            AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(it[j] == g.adjacencyFromVertex(v, j));
            }
        }
        {
            // operator*
            AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(*it == g.adjacencyFromVertex(v, j));
                ++it;
            }
            test(it == g.adjacenciesFromVertexEnd(v));
        }
        // adjacenciesToVertex
        {
            // operator[]
            AdjacencyIterator it = g.adjacenciesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(it[j] == g.adjacencyToVertex(v, j));
            }
        }
        {
            // operator*, operator++, operator==
            AdjacencyIterator it = g.adjacenciesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
                test(*it == g.adjacencyToVertex(v, j));
                ++it;
            }
            test(it == g.adjacenciesToVertexEnd(v));
        }
    }

    // operator==, operator!=
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
        for(std::size_t j = 0; j < g.numberOfVertices() - 1; ++j) {
            for(std::size_t w = 0; w < g.numberOfVertices(); ++w) {
                AdjacencyIterator it2 = g.adjacenciesFromVertexBegin(w);
                for(std::size_t k = 0; k < g.numberOfVertices() - 1; ++k) {
                    test((it + j == it2 + k) == (v == w && j == k));
                    test((it + j != it2 + k) == (v != w || j != k));
                }
            }
        }
    }

    // increment and decrement operators (consistency with operator==, operator!=)
    // TODO: test edge indices
    {
        AdjacencyIterator it = g.adjacenciesFromVertexBegin(1);
        // operator++ (prefix)
        AdjacencyIterator it2 = ++it;
        test(it2 == it);
        test((*it).vertex() == 2);
        test((*it).edge() == g.findEdge(1, 2).second);
        // operator-- (prefix)
        it2 = --it;
        test(it2 == it);
        test((*it).vertex() == 0);
        test((*it).edge() == g.findEdge(1, 0).second);
        // operator++ (postfix)
        it2 = it++;
        test(it2 != it);
        test((*it).vertex() == 2);
        test((*it).edge() == g.findEdge(1, 2).second);
        test((*it2).vertex() == 0);
        test((*it2).edge() == g.findEdge(1, 0).second);
        // operator-- (postfix)
        it2 = it--;
        test(it2 != it);
        test((*it).vertex() == 0);
        test((*it).edge() == g.findEdge(1, 0).second);
        test((*it2).vertex() == 2);
        test((*it2).edge() == g.findEdge(1, 2).second);
        // operator+=
        it2 = (it += 2);
        test(it2 == it);
        test((*it).vertex() == 3);
        test((*it).edge() == g.findEdge(1, 3).second);
        // operator-=
        it2 = (it -= 2);
        test(it2 == it);
        test((*it).vertex() == 0);
        test((*it).edge() == g.findEdge(1, 0).second);
    }

    // operator+, operator-
    // TODO: test edge indices
    {
        AdjacencyIterator it = g.adjacenciesFromVertexBegin(1);
        // operator+
        AdjacencyIterator it2 = it + 2;
        test(it2 != it);
        test((*it).vertex() == 0);
        test((*it).edge() == g.findEdge(1, 0).second);
        test((*it2).vertex() == 3);
        test((*it2).edge() == g.findEdge(1, 3).second);
        // operator-
        AdjacencyIterator it3 = it2 - 2;
        test(it3 != it2);
        test((*it2).vertex() == 3);
        test((*it2).edge() == g.findEdge(1, 3).second);
        test((*it3).vertex() == 0);
        test((*it3).edge() == g.findEdge(1, 0).second);
    }

    // operator>, operator>=, operator<, operator<=
    {
        AdjacencyIterator it = g.adjacenciesFromVertexBegin(1);
        AdjacencyIterator it2 = it;
        ++it2;
        test(!(it2 < it));
        test(!(it2 <= it));
        test(it2 > it);
        test(it2 >= it);
    }
}

int main() {
    testConstructionAndNumbers(); // explicit test
    testFindEdge(); // explicit test
    testAdjacency();
    testVertexIterator();
    testEdgeIterator();
    testAdjacencyIterator();

    return 0;
}
