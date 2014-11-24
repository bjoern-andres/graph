#include <stdexcept>

#include "andres/graph/graph.hxx"

inline void test(const bool& pred) { 
    if(!pred) throw std::runtime_error("Test failed."); 
}

void testConstruction() {
    {
        andres::graph::Graph<> g;
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }
    {
        size_t n = 10;
        andres::graph::Graph<> g(n);
        test(!g.multipleEdgesEnabled());
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

void testAssign() {
    {
        andres::graph::Graph<> g;
        g.reserveVertices(5);
        g.reserveEdges(10);
        g.multipleEdgesEnabled() = true;
        g.assign();
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }
    {
        andres::graph::Graph<> g;
        g.reserveVertices(5);
        g.reserveEdges(10);
        g.multipleEdgesEnabled() = true;
        g.assign(10);
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 10);
        test(g.numberOfEdges() == 0);
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 10);
        test(g.numberOfEdges() == 0);
    }
    {
        size_t n = 10;
        andres::graph::Graph<> g(n);
        g.reserveVertices(5);
        g.reserveEdges(10);
        g.multipleEdgesEnabled() = true;
        g.assign();
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }
    {
        size_t n = 10;
        andres::graph::Graph<> g(n);
        g.reserveVertices(5);
        g.reserveEdges(10);
        g.multipleEdgesEnabled() = true;
        g.assign(n);
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 10);
        test(g.numberOfEdges() == 0);
        g.reserveVertices(5);
        g.reserveEdges(10);
        test(g.numberOfVertices() == 10);
        test(g.numberOfEdges() == 0);
    }
}

void testVertexInsertion() {
    andres::graph::Graph<> g;

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
    andres::graph::Graph<> g(3);
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

    test(g.numberOfVertices() == 3);
    test(g.numberOfEdges() == 9);
    test(g.vertexOfEdge(0, 0) == 1);
    test(g.vertexOfEdge(0, 1) == 2);
    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 1);
    test(g.vertexOfEdge(2, 0) == 1);
    test(g.vertexOfEdge(2, 1) == 1);
    test(g.vertexOfEdge(3, 0) == 1);
    test(g.vertexOfEdge(3, 1) == 1);
    test(g.vertexOfEdge(4, 0) == 0);
    test(g.vertexOfEdge(4, 1) == 2);
    test(g.vertexOfEdge(5, 0) == 0);
    test(g.vertexOfEdge(5, 1) == 2);
    test(g.vertexOfEdge(6, 0) == 2);
    test(g.vertexOfEdge(6, 1) == 2);
    test(g.vertexOfEdge(7, 0) == 2);
    test(g.vertexOfEdge(7, 1) == 2);
    test(g.vertexOfEdge(8, 0) == 0);
    test(g.vertexOfEdge(8, 1) == 1);

    // test adjacencies of vertex 0
    {
        test(g.numberOfEdgesFromVertex(0) == 4);
        test(g.numberOfEdgesToVertex(0) == 4);
        
        test(g.vertexFromVertex(0, 0) == 1);        
        test(g.vertexToVertex(0, 0) == 1);              
        test(g.edgeFromVertex(0, 0) == 1);
        test(g.edgeToVertex(0, 0) == 1);
        test(g.adjacencyFromVertex(0, 0).vertex() == 1);
        test(g.adjacencyFromVertex(0, 0).edge() == 1);
        test(g.adjacencyToVertex(0, 0).vertex() == 1);
        test(g.adjacencyToVertex(0, 0).edge() == 1);

        test(g.vertexFromVertex(0, 1) == 1);
        test(g.vertexToVertex(0, 1) == 1);        
        test(g.edgeFromVertex(0, 1) == 8);
        test(g.edgeToVertex(0, 1) == 8);
        test(g.adjacencyFromVertex(0, 1).vertex() == 1);
        test(g.adjacencyFromVertex(0, 1).edge() == 8);
        test(g.adjacencyToVertex(0, 1).vertex() == 1);
        test(g.adjacencyToVertex(0, 1).edge() == 8);

        test(g.vertexFromVertex(0, 2) == 2);
        test(g.vertexToVertex(0, 2) == 2);        
        test(g.edgeFromVertex(0, 2) == 4);
        test(g.edgeToVertex(0, 2) == 4);
        test(g.adjacencyFromVertex(0, 2).vertex() == 2);
        test(g.adjacencyFromVertex(0, 2).edge() == 4);
        test(g.adjacencyToVertex(0, 2).vertex() == 2);
        test(g.adjacencyToVertex(0, 2).edge() == 4);

        test(g.vertexFromVertex(0, 3) == 2);
        test(g.vertexToVertex(0, 3) == 2);        
        test(g.edgeFromVertex(0, 3) == 5);
        test(g.edgeToVertex(0, 3) == 5);
        test(g.adjacencyFromVertex(0, 3).vertex() == 2);
        test(g.adjacencyFromVertex(0, 3).edge() == 5);
        test(g.adjacencyToVertex(0, 3).vertex() == 2);
        test(g.adjacencyToVertex(0, 3).edge() == 5);
    }

    // test adjacencies of vertex 1
    {
        test(g.numberOfEdgesFromVertex(1) == 5);
        test(g.numberOfEdgesToVertex(1) == 5);

        test(g.vertexFromVertex(1, 0) == 0);        
        test(g.vertexToVertex(1, 0) == 0);        
        test(g.edgeFromVertex(1, 0) == 1);
        test(g.edgeToVertex(1, 0) == 1);
        test(g.adjacencyFromVertex(1, 0).vertex() == 0);
        test(g.adjacencyFromVertex(1, 0).edge() == 1);
        test(g.adjacencyToVertex(1, 0).vertex() == 0);
        test(g.adjacencyToVertex(1, 0).edge() == 1);

        test(g.vertexFromVertex(1, 1) == 0);
        test(g.vertexToVertex(1, 1) == 0);
        test(g.edgeFromVertex(1, 1) == 8);
        test(g.edgeToVertex(1, 1) == 8);
        test(g.adjacencyFromVertex(1, 1).vertex() == 0);
        test(g.adjacencyFromVertex(1, 1).edge() == 8);
        test(g.adjacencyToVertex(1, 1).vertex() == 0);
        test(g.adjacencyToVertex(1, 1).edge() == 8);

        test(g.vertexFromVertex(1, 2) == 1);
        test(g.vertexToVertex(1, 2) == 1);
        test(g.edgeFromVertex(1, 2) == 2);
        test(g.edgeToVertex(1, 2) == 2);
        test(g.adjacencyFromVertex(1, 2).vertex() == 1);
        test(g.adjacencyFromVertex(1, 2).edge() == 2);
        test(g.adjacencyToVertex(1, 2).vertex() == 1);
        test(g.adjacencyToVertex(1, 2).edge() == 2);

        test(g.vertexFromVertex(1, 3) == 1);
        test(g.vertexToVertex(1, 3) == 1);
        test(g.edgeFromVertex(1, 3) == 3);
        test(g.edgeToVertex(1, 3) == 3);
        test(g.adjacencyFromVertex(1, 3).vertex() == 1);
        test(g.adjacencyFromVertex(1, 3).edge() == 3);
        test(g.adjacencyToVertex(1, 3).vertex() == 1);
        test(g.adjacencyToVertex(1, 3).edge() == 3);

        test(g.vertexFromVertex(1, 4) == 2);
        test(g.vertexToVertex(1, 4) == 2);
        test(g.edgeFromVertex(1, 4) == 0);
        test(g.edgeToVertex(1, 4) == 0);
        test(g.adjacencyFromVertex(1, 4).vertex() == 2);
        test(g.adjacencyFromVertex(1, 4).edge() == 0);
        test(g.adjacencyToVertex(1, 4).vertex() == 2);
        test(g.adjacencyToVertex(1, 4).edge() == 0);
    }

    // test adjacencies of vertex 2
    {
        test(g.numberOfEdgesFromVertex(2) == 5);
        test(g.numberOfEdgesToVertex(2) == 5);

        test(g.vertexFromVertex(2, 0) == 0);        
        test(g.vertexToVertex(2, 0) == 0);        
        test(g.edgeFromVertex(2, 0) == 4);
        test(g.edgeToVertex(2, 0) == 4);
        test(g.adjacencyFromVertex(2, 0).vertex() == 0);
        test(g.adjacencyFromVertex(2, 0).edge() == 4);
        test(g.adjacencyToVertex(2, 0).vertex() == 0);
        test(g.adjacencyToVertex(2, 0).edge() == 4);
        
        test(g.vertexFromVertex(2, 1) == 0);
        test(g.vertexToVertex(2, 1) == 0);
        test(g.edgeFromVertex(2, 1) == 5);
        test(g.edgeToVertex(2, 1) == 5);
        test(g.adjacencyFromVertex(2, 1).vertex() == 0);
        test(g.adjacencyFromVertex(2, 1).edge() == 5);
        test(g.adjacencyToVertex(2, 1).vertex() == 0);
        test(g.adjacencyToVertex(2, 1).edge() == 5);
        
        test(g.vertexFromVertex(2, 2) == 1);
        test(g.vertexToVertex(2, 2) == 1);
        test(g.edgeFromVertex(2, 2) == 0);
        test(g.edgeToVertex(2, 2) == 0);
        test(g.adjacencyFromVertex(2, 2).vertex() == 1);
        test(g.adjacencyFromVertex(2, 2).edge() == 0);
        test(g.adjacencyToVertex(2, 2).vertex() == 1);
        test(g.adjacencyToVertex(2, 2).edge() == 0);
        
        test(g.vertexFromVertex(2, 3) == 2);
        test(g.vertexToVertex(2, 3) == 2);
        test(g.edgeFromVertex(2, 3) == 6);
        test(g.edgeToVertex(2, 3) == 6);
        test(g.adjacencyFromVertex(2, 3).vertex() == 2);
        test(g.adjacencyFromVertex(2, 3).edge() == 6);
        test(g.adjacencyToVertex(2, 3).vertex() == 2);
        test(g.adjacencyToVertex(2, 3).edge() == 6);
        
        test(g.vertexFromVertex(2, 4) == 2);
        test(g.vertexToVertex(2, 4) == 2);
        test(g.edgeFromVertex(2, 4) == 7);
        test(g.edgeToVertex(2, 4) == 7);
        test(g.adjacencyFromVertex(2, 4).vertex() == 2);
        test(g.adjacencyFromVertex(2, 4).edge() == 7);
        test(g.adjacencyToVertex(2, 4).vertex() == 2);
        test(g.adjacencyToVertex(2, 4).edge() == 7);
    }
}

void testIterators() {
    andres::graph::Graph<> g(4);
    g.insertEdge(0, 1); // 0
    g.insertEdge(1, 2); // 1
    g.insertEdge(2, 3); // 2
    g.insertEdge(3, 0); // 3

    // vertex iterator, from vertex
    {
        andres::graph::Graph<>::VertexIterator it = g.verticesFromVertexBegin(0);
        test(*it == 1);
        test(it[0] == 1);
        test(it[1] == 3);
        ++it;
        test(*it == 3);
        test(it[0] == 3);
        ++it;
        test(it == g.verticesFromVertexEnd(0));
        test(it-- == g.verticesFromVertexEnd(0));
        test(--it == g.verticesFromVertexBegin(0));
        test((it += 2) == g.verticesFromVertexEnd(0));
        test((it -= 2) == g.verticesFromVertexBegin(0));
        test(it++ == g.verticesFromVertexBegin(0));
        test(++it == g.verticesFromVertexEnd(0));
        test(it - 2 == g.verticesFromVertexBegin(0));
        it -= 2;
        test(it + 2 == g.verticesFromVertexEnd(0));
    }

    // vertex iterator, to vertex
    {
        andres::graph::Graph<>::VertexIterator it = g.verticesToVertexBegin(0);
        test(*it == 1);
        test(it[0] == 1);
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
        andres::graph::Graph<>::EdgeIterator it = g.edgesFromVertexBegin(0);
        test(*it == 0);
        test(it[0] == 0);
        test(it[1] == 3);
        ++it;
        test(*it == 3);
        test(it[0] == 3);
        ++it;
        test(it == g.edgesFromVertexEnd(0));
        test(it-- == g.edgesFromVertexEnd(0));
        test(--it == g.edgesFromVertexBegin(0));
        test((it += 2) == g.edgesFromVertexEnd(0));
        test((it -= 2) == g.edgesFromVertexBegin(0));
        test(it++ == g.edgesFromVertexBegin(0));
        test(++it == g.edgesFromVertexEnd(0));
        test(it - 2 == g.edgesFromVertexBegin(0));
        it -= 2;
        test(it + 2 == g.edgesFromVertexEnd(0));
    }

    // edge iterator, to vertex
    {
        andres::graph::Graph<>::EdgeIterator it = g.edgesToVertexBegin(0);
        test(*it == 0);
        test(it[0] == 0);
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
        andres::graph::Graph<>::AdjacencyIterator it = g.adjacenciesFromVertexBegin(0);
        test(it->vertex() == 1);
        test(it[0].vertex() == 1);
        test(it[1].vertex() == 3);
        ++it;
        test(it->vertex() == 3);
        test(it[0].vertex() == 3);
        ++it;
        test(it == g.verticesFromVertexEnd(0));
        test(it-- == g.verticesFromVertexEnd(0));
        test(--it == g.verticesFromVertexBegin(0));
        test((it += 2) == g.verticesFromVertexEnd(0));
        test((it -= 2) == g.verticesFromVertexBegin(0));
        test(it++ == g.verticesFromVertexBegin(0));
        test(++it == g.verticesFromVertexEnd(0));
        test(it - 2 == g.verticesFromVertexBegin(0));
        it -= 2;
        test(it + 2 == g.verticesFromVertexEnd(0));
    }

    // adjacency iterator, to vertex
    {
        andres::graph::Graph<>::AdjacencyIterator it = g.adjacenciesToVertexBegin(0);
        test(it->vertex() == 1);
        test(it[0].vertex() == 1);
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

template<typename I>
class IteratorTest {
public:
    typedef I iterator;
    typedef typename std::iterator_traits<iterator>::value_type value_type;
    typedef typename std::iterator_traits<iterator>::difference_type difference_type;
    typedef typename std::iterator_traits<iterator>::pointer pointer;
    typedef typename std::iterator_traits<iterator>::reference reference;
    typedef typename std::iterator_traits<iterator>::iterator_category iterator_category;
    
    void test(iterator begin,iterator end) {
        difference_type d = std::distance(begin,end);
        assert(d>0);
        iterator it0 = std::next(begin);
        {
            iterator it1 = begin;
            std::advance(it1,d);
            ::test(it1!=begin);
            ::test(it1==end);
        }
        {
            // Decrement it2 by d-1
            iterator it1 = std::prev(end,d-1);
            ::test(it1!=end);
            ::test(it1==it0);
        }
    }
};

void testIteratorCompiling() {
    typedef andres::graph::Graph<> Graph;
    typedef typename Graph::AdjacencyIterator AdjacencyIterator ;
    typedef typename Graph::VertexIterator VertexIterator;
    typedef typename Graph::EdgeIterator EdgeIterator;
    
    Graph g(4);
    g.insertEdge(0, 1); // 0
    g.insertEdge(1, 2); // 1
    g.insertEdge(2, 3); // 2
    g.insertEdge(3, 0); // 3

    std::size_t pivot = 1;
    
    {
        VertexIterator begin = g.verticesFromVertexBegin(pivot);
        VertexIterator end = g.verticesFromVertexEnd(pivot);
        IteratorTest<VertexIterator> vertexIteratorTest;
        vertexIteratorTest.test(begin,end);
    }
    {
        EdgeIterator begin = g.edgesFromVertexBegin(pivot);
        EdgeIterator end = g.edgesFromVertexEnd(pivot);
        IteratorTest<EdgeIterator> edgeIteratorTest;
        edgeIteratorTest.test(begin,end);
    }
    {
        AdjacencyIterator begin = g.adjacenciesFromVertexBegin(pivot);
        AdjacencyIterator end = g.adjacenciesFromVertexEnd(pivot);
        IteratorTest<AdjacencyIterator> adjacencyIteratorTest;
        adjacencyIteratorTest.test(begin,end);
    }
}

void testEdgeRemoval() {
    andres::graph::Graph<> g(3);
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

    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 2);

    test(g.vertexOfEdge(2, 0) == 1);
    test(g.vertexOfEdge(2, 1) == 1);
}

void testVertexRemoval() 
{
    andres::graph::Graph<> g(3);
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

    test(g.vertexOfEdge(0, 0) == 0);
    test(g.vertexOfEdge(0, 1) == 1);

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
    andres::graph::Graph<> g(4);
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
    test(p.first == true);
    test(p.second == 3);

    p = g.findEdge(1, 0);
    test(p.first == true);
    test(p.second == 0);

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
    test(p.first == true);
    test(p.second == 1);

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
    test(p.first == true);
    test(p.second == 2);

    p = g.findEdge(3, 3);
    test(p.first == false);
}

void testMultipleEdges() {
    andres::graph::Graph<> g(4);
    g.insertEdge(0, 1); 
    test(g.numberOfEdges() == 1);
    test(g.insertEdge(0, 1) == 0);
    test(g.numberOfEdges() == 1);

    g.multipleEdgesEnabled() = true;
    g.insertEdge(0, 1); 
    test(g.numberOfEdges() == 2);
    test(g.vertexOfEdge(0, 0) == 0);
    test(g.vertexOfEdge(0, 1) == 1);
    test(g.vertexOfEdge(1, 0) == 0);
    test(g.vertexOfEdge(1, 1) == 1);
}

int main() {
    testConstruction();
    testAssign();
    testVertexInsertion();
    testIterators();
    testIteratorCompiling();
    testEdgeInsertion();
    testEdgeRemoval();
    testVertexRemoval();
    testfindEdge();
    testMultipleEdges();

    return 0;
}
