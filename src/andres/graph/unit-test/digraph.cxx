#include <stdexcept>

#include "andres/graph/digraph.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

void testConstruction() {
    {
        andres::graph::Digraph<> g;
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }
    {
        size_t n = 10;
        andres::graph::Digraph<> g(n);
        test(g.numberOfVertices() == 10);
        test(g.numberOfEdges() == 0);
        for(size_t j=0; j<n; ++j) {
            test(g.numberOfEdgesFromVertex(j) == 0);
            test(g.numberOfEdgesToVertex(j) == 0);
        }
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 10);
        test(g.numberOfEdges() == 0);
    }
}

void testVertexInsertion() {
    andres::graph::Digraph<> g;

    size_t index = g.insertVertex();
    test(index == 0);
    test(g.numberOfVertices() == 1);
    test(g.numberOfEdges() == 0);
    test(g.numberOfEdgesFromVertex(0) == 0);
    test(g.numberOfEdgesToVertex(0) == 0);

    index = g.insertVertex();
    test(index == 1);
    test(g.numberOfVertices() == 2);
    test(g.numberOfEdges() == 0);
    for(size_t j=0; j<2; ++j) {
        test(g.numberOfEdgesFromVertex(j) == 0);
        test(g.numberOfEdgesToVertex(j) == 0);
    }

    index = g.insertVertices(5);
    test(index == 2);
    test(g.numberOfVertices() == 7);
    test(g.numberOfEdges() == 0);
    for(size_t j=0; j<7; ++j) {
        test(g.numberOfEdgesFromVertex(j) == 0);
        test(g.numberOfEdgesToVertex(j) == 0);
    }
}

void testEdgeInsertion() {
    andres::graph::Digraph<> g(4);
    g.insertEdge(0, 1); // 0
    g.insertEdge(1, 2); // 1
    g.insertEdge(2, 3); // 2
    g.insertEdge(3, 0); // 3
    g.insertEdge(2, 0); // 4

    test(g.numberOfVertices() == 4);
    test(g.numberOfEdges() == 5);

    test(g.vertexOfEdge(0, 0) == 0);
    test(g.vertexOfEdge(0, 1) == 1);

    test(g.vertexOfEdge(1, 0) == 1);
    test(g.vertexOfEdge(1, 1) == 2);

    test(g.vertexOfEdge(2, 0) == 2);
    test(g.vertexOfEdge(2, 1) == 3);

    test(g.vertexOfEdge(3, 0) == 3);
    test(g.vertexOfEdge(3, 1) == 0);

    test(g.vertexOfEdge(4, 0) == 2);
    test(g.vertexOfEdge(4, 1) == 0);

    test(g.numberOfEdgesFromVertex(0) == 1);
    test(g.vertexFromVertex(0, 0) == 1);
    test(g.adjacencyFromVertex(0, 0).vertex() == 1);
    test(g.adjacencyFromVertex(0, 0).edge() == 0);
    test(g.numberOfEdgesToVertex(0) == 2);       
    test(g.vertexToVertex(0, 0) == 2);
    test(g.vertexToVertex(0, 1) == 3);
    test(g.adjacencyToVertex(0, 0).vertex() == 2);
    test(g.adjacencyToVertex(0, 0).edge() == 4);
    test(g.adjacencyToVertex(0, 1).vertex() == 3);
    test(g.adjacencyToVertex(0, 1).edge() == 3);
        
    test(g.numberOfEdgesFromVertex(1) == 1);
    test(g.vertexFromVertex(1, 0) == 2);
    test(g.adjacencyFromVertex(1, 0).vertex() == 2);
    test(g.adjacencyFromVertex(1, 0).edge() == 1);
    test(g.numberOfEdgesToVertex(1) == 1);
    test(g.vertexToVertex(1, 0) == 0);
    test(g.adjacencyToVertex(1, 0).vertex() == 0);
    test(g.adjacencyToVertex(1, 0).edge() == 0);

    test(g.numberOfEdgesFromVertex(2) == 2);
    test(g.vertexFromVertex(2, 0) == 0);
    test(g.vertexFromVertex(2, 1) == 3);
    test(g.adjacencyFromVertex(2, 0).vertex() == 0);
    test(g.adjacencyFromVertex(2, 0).edge() == 4);
    test(g.adjacencyFromVertex(2, 1).vertex() == 3);
    test(g.adjacencyFromVertex(2, 1).edge() == 2);
    test(g.numberOfEdgesToVertex(2) == 1);
    test(g.vertexToVertex(2, 0) == 1);
    test(g.adjacencyToVertex(2, 0).vertex() == 1);
    test(g.adjacencyToVertex(2, 0).edge() == 1);

    test(g.numberOfEdgesFromVertex(3) == 1);
    test(g.vertexFromVertex(3, 0) == 0);
    test(g.adjacencyFromVertex(3, 0).vertex() == 0);
    test(g.adjacencyFromVertex(3, 0).edge() == 3);
    test(g.numberOfEdgesToVertex(3) == 1);
    test(g.vertexToVertex(3, 0) == 2);
    test(g.adjacencyToVertex(3, 0).vertex() == 2);
    test(g.adjacencyToVertex(3, 0).edge() == 2);
}

void testIterators() {
    andres::graph::Digraph<> g(4);
    g.insertEdge(0, 1); // 0
    g.insertEdge(1, 2); // 1
    g.insertEdge(2, 3); // 2
    g.insertEdge(3, 0); // 3
    g.insertEdge(2, 0); // 4

    // vertex iterator, from vertex
    {
        andres::graph::Digraph<>::VertexIterator it = g.verticesFromVertexBegin(2);
        test(*it == 0);
        test(it[0] == 0);
        test(it[1] == 3);
        ++it;
        test(*it == 3);
        test(it[0] == 3);
        ++it;
        test(it == g.verticesFromVertexEnd(2));
        test(it-- == g.verticesFromVertexEnd(2));
        test(--it == g.verticesFromVertexBegin(2));
        test((it += 2) == g.verticesFromVertexEnd(2));
        test((it -= 2) == g.verticesFromVertexBegin(2));
        test(it++ == g.verticesFromVertexBegin(2));
        test(++it == g.verticesFromVertexEnd(2));
        test(it - 2 == g.verticesFromVertexBegin(2));
        it -= 2;
        test(it + 2 == g.verticesFromVertexEnd(2));
    }

    // vertex iterator, to vertex
    {
        andres::graph::Digraph<>::VertexIterator it = g.verticesToVertexBegin(0);
        test(*it == 2);
        test(it[0] == 2);
        test(it[1] == 3);
        ++it;
        test(*it == 3);
        test(it[0] == 3);
        ++it;
        test(it == g.verticesToVertexEnd(0));
        test(it-- == g.verticesToVertexEnd(0));
        test(--it == g.verticesToVertexBegin(0));
        test((it += 2) == g.verticesToVertexEnd(0));
        test((it -= 2) == g.verticesToVertexBegin(0));
        test(it++ == g.verticesToVertexBegin(0));
        test(++it == g.verticesToVertexEnd(0));
        test(it - 2 == g.verticesToVertexBegin(0));
        it -= 2;
        test(it + 2 == g.verticesToVertexEnd(0));
    }

    // edge iterator, from vertex
    {
        andres::graph::Digraph<>::EdgeIterator it = g.edgesFromVertexBegin(2);
        std::cout << *it << std::endl;
        test(*it == 4);
        test(it[0] == 4);
        test(it[1] == 2);
        ++it;
        test(*it == 2);
        test(it[0] == 2);
        ++it;
        test(it == g.edgesFromVertexEnd(2));
        test(it-- == g.edgesFromVertexEnd(2));
        test(--it == g.edgesFromVertexBegin(2));
        test((it += 2) == g.edgesFromVertexEnd(2));
        test((it -= 2) == g.edgesFromVertexBegin(2));
        test(it++ == g.edgesFromVertexBegin(2));
        test(++it == g.edgesFromVertexEnd(2));
        test(it - 2 == g.edgesFromVertexBegin(2));
        it -= 2;
        test(it + 2 == g.edgesFromVertexEnd(2));
    }

    // edge iterator, to vertex
    {
        andres::graph::Digraph<>::EdgeIterator it = g.edgesToVertexBegin(0);
        test(*it == 4);
        test(it[0] == 4);
        test(it[1] == 3);
        ++it;
        test(*it == 3);
        test(it[0] == 3);
        ++it;
        test(it == g.edgesToVertexEnd(0));
        test(it-- == g.edgesToVertexEnd(0));
        test(--it == g.edgesToVertexBegin(0));
        test((it += 2) == g.edgesToVertexEnd(0));
        test((it -= 2) == g.edgesToVertexBegin(0));
        test(it++ == g.edgesToVertexBegin(0));
        test(++it == g.edgesToVertexEnd(0));
        test(it - 2 == g.edgesToVertexBegin(0));
        it -= 2;
        test(it + 2 == g.edgesToVertexEnd(0));
    }

    // adjacency iterator, from vertex
    {
        andres::graph::Digraph<>::AdjacencyIterator it = 
            g.adjacenciesFromVertexBegin(2);
        test(it->vertex() == 0);
        test(it[0].vertex() == 0);
        test(it[1].vertex() == 3);
        ++it;
        test(it->vertex() == 3);
        test(it[0].vertex() == 3);
        ++it;
        test(it == g.verticesFromVertexEnd(2));
        test(it-- == g.verticesFromVertexEnd(2));
        test(--it == g.verticesFromVertexBegin(2));
        test((it += 2) == g.verticesFromVertexEnd(2));
        test((it -= 2) == g.verticesFromVertexBegin(2));
        test(it++ == g.verticesFromVertexBegin(2));
        test(++it == g.verticesFromVertexEnd(2));
        test(it - 2 == g.verticesFromVertexBegin(2));
        it -= 2;
        test(it + 2 == g.verticesFromVertexEnd(2));
    }

    // adjacency iterator, to vertex
    {
        andres::graph::Digraph<>::AdjacencyIterator it = 
            g.adjacenciesToVertexBegin(0);
        test(it->vertex() == 2);
        test(it[0].vertex() == 2);
        test(it[1].vertex() == 3);
        ++it;
        test(it->vertex() == 3);
        test(it[0].vertex() == 3);
        ++it;
        test(it == g.verticesToVertexEnd(0));
        test(it-- == g.verticesToVertexEnd(0));
        test(--it == g.verticesToVertexBegin(0));
        test((it += 2) == g.verticesToVertexEnd(0));
        test((it -= 2) == g.verticesToVertexBegin(0));
        test(it++ == g.verticesToVertexBegin(0));
        test(++it == g.verticesToVertexEnd(0));
        test(it - 2 == g.verticesToVertexBegin(0));
        it -= 2;
        test(it + 2 == g.verticesToVertexEnd(0));
    }
}

void testEdgeRemoval() {
    andres::graph::Digraph<> g(3);
    g.multipleEdgesEnabled() = true;
    g.insertEdge(1, 2); // 0
    g.insertEdge(1, 0); // 1
    g.insertEdge(1, 1); // 2
    g.insertEdge(1, 1); // 3
    g.insertEdge(0, 2); // 4
    g.insertEdge(2, 0); // 5
    g.insertEdge(2, 2); // 6
    g.insertEdge(2, 2); // 7
    g.insertEdge(0, 1); // 8

    g.eraseEdge(0); // orig. 0
    g.eraseEdge(4); // orig. 4
    g.eraseEdge(4); // orig. 7
    g.eraseEdge(3); // orig. 3
    g.eraseEdge(0); // orig. 8
    g.eraseEdge(1); // orig. 1

    test(g.numberOfVertices() == 3);
    test(g.numberOfEdges() == 3);

    test(g.vertexOfEdge(0, 0) == 2);
    test(g.vertexOfEdge(0, 1) == 2);

    test(g.vertexOfEdge(1, 0) == 2);
    test(g.vertexOfEdge(1, 1) == 0);

    test(g.vertexOfEdge(2, 0) == 1);
    test(g.vertexOfEdge(2, 1) == 1);
}

void testVertexRemoval() 
{
    andres::graph::Digraph<> g(3);
    g.multipleEdgesEnabled() = true;
    g.insertEdge(1, 2); // 0
    g.insertEdge(1, 0); // 1
    g.insertEdge(1, 1); // 2
    g.insertEdge(1, 1); // 3
    g.insertEdge(0, 2); // 4
    g.insertEdge(2, 0); // 5
    g.insertEdge(2, 2); // 6
    g.insertEdge(2, 2); // 7
    g.insertEdge(0, 1); // 8

    // first removal
    g.eraseVertex(0); // orig. 0

    test(g.numberOfVertices() == 2);
    test(g.numberOfEdges() == 5);

    test(g.vertexOfEdge(0, 0) == 1);
    test(g.vertexOfEdge(0, 1) == 0);

    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 0);

    test(g.vertexOfEdge(2, 0) == 1);
    test(g.vertexOfEdge(2, 1) == 1);

    test(g.vertexOfEdge(3, 0) == 1);
    test(g.vertexOfEdge(3, 1) == 1);

    test(g.vertexOfEdge(4, 0) == 0);
    test(g.vertexOfEdge(4, 1) == 0);

    // second removal
    g.eraseVertex(0); // orig. 2

    test(g.numberOfVertices() == 1);
    test(g.numberOfEdges() == 2);

    test(g.vertexOfEdge(0, 0) == 0);
    test(g.vertexOfEdge(0, 1) == 0);

    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 0);
}

void testfindEdge() {
    andres::graph::Digraph<> g(4);
    g.insertEdge(0, 1); // 0
    g.insertEdge(1, 2); // 1
    g.insertEdge(2, 3); // 2
    g.insertEdge(3, 0); // 3

    std::pair<bool, size_t> p;
    p = g.findEdge(0, 0);
    test(p.first == false);

    p = g.findEdge(0, 1);
    test(p.first == true);
    test(p.second == 0);

    p = g.findEdge(0, 2);
    test(p.first == false);

    p = g.findEdge(0, 3);
    test(p.first == false);

    p = g.findEdge(1, 0);
    test(p.first == false);

    p = g.findEdge(1, 1);
    test(p.first == false);

    p = g.findEdge(1, 2);
    test(p.first == true);
    test(p.second == 1);

    p = g.findEdge(1, 3);
    test(p.first == false);

    p = g.findEdge(2, 0);
    test(p.first == false);

    p = g.findEdge(2, 1);
    test(p.first == false);

    p = g.findEdge(2, 2);
    test(p.first == false);

    p = g.findEdge(2, 3);
    test(p.first == true);
    test(p.second == 2);

    p = g.findEdge(3, 0);
    test(p.first == true);
    test(p.second == 3);

    p = g.findEdge(3, 1);
    test(p.first == false);

    p = g.findEdge(3, 2);
    test(p.first == false);

    p = g.findEdge(3, 3);
    test(p.first == false);
}

void testMultipleEdges() {
    andres::graph::Digraph<> g(4);
    g.insertEdge(0, 1); 
    test(g.numberOfEdges() == 1);
    test(g.insertEdge(0, 1) == 0);
    test(g.numberOfEdges() == 1);
    
    g.multipleEdgesEnabled() = true;
    test(g.insertEdge(0, 1) == 1); 
    test(g.numberOfEdges() == 2);
    test(g.vertexOfEdge(0, 0) == 0);
    test(g.vertexOfEdge(0, 1) == 1);
    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 1);
    
    g.multipleEdgesEnabled() = false;
    g.insertEdge(1, 0);
    test(g.vertexOfEdge(0, 0) == 0);
    test(g.vertexOfEdge(0, 1) == 1);
    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 1);
    test(g.vertexOfEdge(2, 0) == 1);
    test(g.vertexOfEdge(2, 1) == 0);
}

int main() {
    testConstruction();
    testVertexInsertion();
    testEdgeInsertion();
    testIterators();
    testEdgeRemoval();
    testVertexRemoval();
    testfindEdge();
    testMultipleEdges();

    return 0;
}
