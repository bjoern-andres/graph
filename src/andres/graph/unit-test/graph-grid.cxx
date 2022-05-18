#include <stdexcept>
#include <iostream>
#include <memory>
#include <algorithm>
#include <vector>
#include <iterator>
#include <type_traits>

#include "andres/graph/grid-graph.hxx"
#include "andres/graph/graph.hxx"

inline void test(const bool& pred) {
    if(!pred)
        throw std::runtime_error("Test failed.");
}

template<unsigned char D, class V0, class V1>
void testGraphsEquals(
    const andres::graph::GridGraph<D, V0> & graph0,
    const andres::graph::GridGraph<D, V1> & graph1
) {
    test(graph1.numberOfVertices() == graph0.numberOfVertices());
    test(graph1.numberOfEdges() == graph0.numberOfEdges());
    for(std::size_t j = 0; j < andres::graph::GridGraph<D, V0>::DIMENSION; ++j) {
        test(graph1.shape(j) == graph0.shape(j));
    }
    for(std::size_t j = 0; j < graph1.numberOfVertices(); ++j)
    for(std::size_t k = 0; k < graph1.numberOfVertices(); ++k) {
        test(graph1.findEdge(j, k).first == graph0.findEdge(j, k).first);
    }
}

#define testThrowsException(PRED, EXC) { \
        {\
            bool threw = false; \
            try { \
                PRED; \
            } catch (EXC& e) {\
                threw = true; \
            } catch (...) {\
                throw std::runtime_error("Test threw different exception type.");\
            }\
            if(!threw) \
                throw std::runtime_error("Test did not throw exception");\
        }\
    }
#define testThrows(PRED) { \
        {\
            bool threw = false; \
            try { \
                PRED; \
            } catch (...) {\
                threw = true; \
            }\
            if(!threw) \
                throw std::runtime_error("Test did not throw exception");\
        }\
    }

typedef andres::graph::Graph<> Graph;
typedef std::array<size_t,2> VertexPair;

struct EdgeEndpoint {
    std::vector<size_t> coord0;
    std::vector<size_t> coord1;
    size_t idx0;
    size_t idx1;
};

struct VertexEntry {
    std::vector<size_t> coord;
    size_t index;
};

// Compute the cumulative product of an input vector
inline size_t cumprod(const std::vector<size_t>& shape, std::vector<size_t>& cumprod) {
    if(shape.size()==0)
        return 0;
    cumprod.resize(shape.size()+1);
    size_t prod = shape[0];
    cumprod[0] = 1;
    {
        size_t i;
        for(i=1; i<shape.size(); ++i) {
            cumprod[i] = prod;
            prod *= shape[i];
        }
        cumprod[i] = prod;
    }
    return prod;
}

// Compute the product of an input vector
inline size_t prod(const std::vector<size_t>& shape) {
    if(shape.size()==0)
        return 0;
    size_t prod = shape[0];
    for(size_t i=1; i<shape.size(); ++i) {
        prod *= shape[i];
    }
    return prod;
}

void createCorrectGridND(const std::vector<size_t>& shape, Graph& g) {
    size_t dimensions = shape.size();
    std::vector<size_t> vertexOffsets;
    size_t numVertices = cumprod(shape, vertexOffsets);
    g.assign();
    g.insertVertices(numVertices);

    std::vector<size_t> maxCoord(dimensions);
    maxCoord.push_back(1); // Terminating condition: Overflow of the last dimension counter: coord == [0, ..., 0 1]
    for(size_t d=0; d<dimensions; ++d) {
        size_t pivotIndex=0;
        std::vector<size_t> edgeShape(shape);
        --edgeShape[d];
        {
            std::vector<size_t> coord(dimensions+1);
            while(coord!=maxCoord) {
                const size_t pairIndex = pivotIndex + vertexOffsets[d];
                {
                    EdgeEndpoint endpoint;
                    endpoint.coord0.insert(endpoint.coord0.begin(),coord.begin(),coord.begin()+dimensions);
                    endpoint.coord1.insert(endpoint.coord1.begin(),coord.begin(),coord.begin()+dimensions);
                    ++endpoint.coord1[d];
                    endpoint.idx0 = pivotIndex;
                    endpoint.idx1 = pairIndex;
                    g.insertEdge(pivotIndex,pairIndex); // insert edge along dimension d
                }
                ++pivotIndex;
                ++coord[0];
                for(size_t i=0; i<dimensions && coord[i]==edgeShape[i]; ++i) {
                    coord[i]=0;
                    ++coord[i+1];
                    if(edgeShape[i]!=shape[i]) { // Border case. Skip line (advance by offset)
                        pivotIndex += vertexOffsets[i];
                    }
                }
            }
        }
    }
}

// Create a LUT with edge endpoints fror all edges in the grid graph of the given shape
void computeEdgeEndpointsND(const std::vector<size_t> shape, std::vector<EdgeEndpoint>& endpoints) {
    size_t dimensions = shape.size();
    std::vector<size_t> vertexOffsets;
    endpoints.clear();

    std::vector<size_t> maxCoord(dimensions);
    maxCoord.push_back(1); // Terminating condition: Overflow of the last dimension counter: coord == [0, ..., 0 1]
    for(size_t d=0; d<dimensions; ++d) {
        size_t pivotIndex=0;
        std::vector<size_t> edgeShape(shape);
        --edgeShape[d];
        {
            std::vector<size_t> coord(dimensions+1);
            while(coord!=maxCoord) {
                const size_t pairIndex = pivotIndex + vertexOffsets[d];
                {
                    EdgeEndpoint endpoint;
                    endpoint.coord0.insert(endpoint.coord0.begin(),coord.begin(),coord.begin()+dimensions);
                    endpoint.coord1.insert(endpoint.coord1.begin(),coord.begin(),coord.begin()+dimensions);
                    ++endpoint.coord1[d];
                    endpoint.idx0 = pivotIndex;
                    endpoint.idx1 = pairIndex;
                    endpoints.push_back(endpoint);
                }
                ++pivotIndex;
                ++coord[0];
                for(size_t i=0; i<dimensions && coord[i]==edgeShape[i]; ++i) {
                    coord[i]=0;
                    ++coord[i+1];
                    if(edgeShape[i]!=shape[i]) { // Border case. Skip line (advance by offset)
                        pivotIndex += vertexOffsets[i];
                    }
                }
            }
        }
    }
}

// Iterate over all coordinates and create a LUT of vertex index coordinates
void computeVertexCoordinatesND(const std::vector<size_t> shape, std::vector<std::vector<size_t> >& vertexCoordinates) {
    const size_t dimensions = shape.size();
    const size_t numVertices = prod(shape);
    {
        const size_t numVertices = prod(shape);
        vertexCoordinates.resize(numVertices);
    }
    std::vector<size_t> coord(dimensions);
    for(size_t idx=0; idx < numVertices ; ++idx) {
        vertexCoordinates[idx] = coord;
        ++coord[0]; // Increase Least Significant VertexCoordinate
        for(size_t i=0; i<dimensions-1 && coord[i]==shape[i]; ++i) { // Propagate
            coord[i]=0;
            ++coord[i+1];
        }
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

    static_assert(std::is_same<pointer,value_type*>::value,"Iterator pointer type mismatch.");
    static_assert(std::is_same<reference,value_type&>::value,"Iterator pointer type mismatch.");
    
    void operator()(iterator begin,iterator end) {
        difference_type d = std::distance(begin,end);
        assert(d>0);
        iterator it0 = std::next(begin);
        {
            iterator it1 = begin;
            std::advance(it1,d);
            test(it1!=begin);
            test(it1==end);
        }
        {
            // Decrement it2 by d-1
            iterator it1 = std::prev(end,d-1);
            test(it1!=end);
            test(it1==it0);
        }
    }
};


template<typename G>
void testIteratorCompile(G& g, size_t pivot) {
    typedef G GridGraph;
    typedef typename GridGraph::AdjacencyIterator AdjacencyIterator ;
    typedef typename GridGraph::VertexIterator VertexIterator;
    typedef typename GridGraph::EdgeIterator EdgeIterator;
    
    {
        VertexIterator begin = g.verticesFromVertexBegin(pivot);
        VertexIterator end = g.verticesFromVertexEnd(pivot);
        IteratorTest<VertexIterator> vertexIteratorTest;
        vertexIteratorTest(begin,end);
    }
    {
        EdgeIterator begin = g.edgesFromVertexBegin(pivot);
        EdgeIterator end = g.edgesFromVertexEnd(pivot);
        IteratorTest<EdgeIterator> edgeIteratorTest;
        edgeIteratorTest(begin,end);
    }
    {
        AdjacencyIterator begin = g.adjacenciesFromVertexBegin(pivot);
        AdjacencyIterator end = g.adjacenciesFromVertexEnd(pivot);
        IteratorTest<AdjacencyIterator> adjacencyIteratorTest;
        adjacencyIteratorTest(begin,end);
    }
}

void testConstructionAndNumbersExplicit() {
    {
        typedef andres::graph::GridGraph<4> GridGraph;
        typedef typename GridGraph::VertexCoordinate VertexCoordinate;
        const andres::graph::GridGraph<4> g(VertexCoordinate( {{4,6,7,3}}));
        test(g.numberOfVertices() == 4*6*7*3);
        test(g.numberOfEdges() == 3*6*7*3+4*5*7*3+4*6*6*3+4*6*7*2);
    }
}

template<typename G>
void testConstructionAndNumbersND(const std::vector<size_t>& shape) {
    typedef G GridGraph;
    typedef typename GridGraph::VertexCoordinate VertexCoordinate;
    {
        const GridGraph g;
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }
    {
        VertexCoordinate vcShape;
        std::copy(shape.begin(),shape.end(),vcShape.begin());
        GridGraph g(vcShape);
        test(g.numberOfVertices() == prod(shape));
        size_t numEdges = 0;
        for(size_t i=0; i<g.DIMENSION; ++i) {
            std::vector<size_t> edgeShape = shape;
            --edgeShape[i];
            numEdges+=prod(edgeShape);
        }
        test(g.numberOfEdges() == numEdges);
    }

    {
        std::vector<std::vector<size_t> > vertexCoordinates;
        GridGraph g;
        {
            std::array<size_t,4> shape_array;
            std::copy(shape.begin(),shape.end(),shape_array.begin());
            g.assign(shape_array);
            for(size_t i=0; i<g.DIMENSION; ++i)
                test(g.shape(i) == shape_array[i]);
        }
        computeVertexCoordinatesND(shape, vertexCoordinates);
        for(size_t idx=0; idx<vertexCoordinates.size(); ++idx) {
            std::vector<size_t>& coord = vertexCoordinates[idx];
            test(!g.multipleEdgesEnabled());
            test(g.numberOfVertices() == prod(shape));
            VertexCoordinate shape;
            for(size_t i=0; i<g.DIMENSION; ++i)
                shape[i] = g.shape(i);
            {
                std::size_t edges = 0;
                for(size_t i=0; i<shape.size(); ++i) {
                    if(coord[i]>0) ++edges;
                    if(coord[i]<shape[i]-1) ++edges;
                }
                test(g.numberOfEdgesFromVertex(idx) == edges);
                test(g.numberOfEdgesToVertex(idx) == edges);
                {
                    VertexCoordinate arrCoord;
                    std::copy(coord.begin(),coord.end(),arrCoord.begin());
                    {
                        VertexCoordinate vCoord;
                        g.vertex(idx,vCoord);
                        test(vCoord == arrCoord);
                    }
                    test(g.vertex(arrCoord) == idx);
                }
            }
        }
    }
}




// tests findEdge, Orientations, Edge mappings
template<typename G>
void testOrientationsND(const std::vector<size_t>& shape) {
    typedef G GridGraph;
    typedef typename GridGraph::VertexCoordinate VertexCoordinate;
    typedef typename GridGraph::EdgeCoordinate EdgeCoordinate;

    GridGraph g;
    {
        VertexCoordinate vcShape;
        std::copy(shape.begin(),shape.end(),vcShape.begin());
        g.assign(vcShape);
    }

    for(size_t u=0; u<g.numberOfVertices(); ++u) {
        VertexCoordinate coordinate;
        g.vertex(u,coordinate);
        for(size_t d=0; d<g.DIMENSION; ++d) {
            if(coordinate[d]>0) {
                const EdgeCoordinate ecSmaller(coordinate,d,true);
                test(coordinate[d]==ecSmaller.pivotCoordinate_[d]+1); // internals
                const size_t edge = g.edge(ecSmaller); // Get Edge index
                size_t target = g.vertexOfEdge(edge,0);
                if(target==u) // request oposite endpoint
                    target = g.vertexOfEdge(edge,1);
                std::pair<bool,size_t> qr = g.findEdge(u,target);
                test(qr.first == true);
                test(qr.second == edge);
                VertexCoordinate vcTarget;
                g.vertex(target,vcTarget);
                test(vcTarget[d]+1==coordinate[d]);
            }
            if(coordinate[d]<shape[d]-1) {
                const EdgeCoordinate ecSmaller(coordinate,d,false);
                test(coordinate[d]==ecSmaller.pivotCoordinate_[d]); // internals
                const size_t edge = g.edge(ecSmaller); // Get Edge index
                size_t target = g.vertexOfEdge(edge,0);
                if(target==u) // request oposite endpoint
                    target = g.vertexOfEdge(edge,1);
                std::pair<bool,size_t> qr = g.findEdge(u,target);
                test(qr.first == true);
                test(qr.second == edge);
                VertexCoordinate vcTarget;
                g.vertex(target,vcTarget);
                test(vcTarget[d]==coordinate[d]+1);
            }
        }
    }
}


// tests consistency of adjacency functions with findEdge
template<typename G>
void testFindEdgeND(const std::vector<size_t>& shape) {
    typedef G GridGraph;
    typedef typename GridGraph::VertexCoordinate VertexCoordinate;
    typedef typename GridGraph::AdjacencyType AdjacencyType;

    GridGraph gg;
    {
        VertexCoordinate vcShape;
        std::copy(shape.begin(),shape.end(),vcShape.begin());
        gg.assign(vcShape);
    }
    Graph g;
    createCorrectGridND(shape,g);

    const size_t INVALID_EDGE = g.numberOfEdges();
    std::vector<size_t> adjacentEdges(g.numberOfVertices()); // This will be valid if the vertex is adjacent
    for(std::size_t s = 0; s < g.numberOfVertices(); ++s) {
        // Set the correct adjacencies
        {
            std::fill(adjacentEdges.begin(),adjacentEdges.end(),INVALID_EDGE);
            size_t numAdjacentVertices = g.numberOfEdgesFromVertex(s);
            for(size_t i=0; i<numAdjacentVertices; ++i) {
                AdjacencyType adjacency = g.adjacencyFromVertex(s,i);
                adjacentEdges[adjacency.vertex()] = adjacency.edge();
            }
        }

        for(std::size_t t = 0; t < g.numberOfVertices(); ++t) {
            std::pair<bool,size_t> queryResult = gg.findEdge(s,t);
            const size_t adjacentEdge = adjacentEdges[t];
            if(adjacentEdge != INVALID_EDGE) {
                test(queryResult.first == true); // Test detection
                test(queryResult.second == adjacentEdge); // test correct recovery
                test(gg.findEdge(t,s) == queryResult); // Test Symmetry
            } else {
                test(queryResult.first == false); // Test detection
                test(gg.findEdge(t,s).first == false); // test symmetry
            }
            // Test insertEdge
            {
                std::pair<bool, std::size_t> qr = gg.findEdge(s,t);
                if(qr.first == true) {
                    test(qr.second == gg.insertEdge(s,t));
                } else {
                    testThrowsException(gg.insertEdge(s,t),std::runtime_error);
                }
            }
        }
    }
}

// tests consistency of adjacency functions with findEdge
template<typename G>
void testAdjacencyND(std::vector<size_t> shape) {
    typedef G GridGraph;
    typedef typename GridGraph::AdjacencyType AdjacencyType;

    GridGraph gg;
    {
        std::array<size_t,4> shape_array;
        std::copy(shape.begin(),shape.end(),shape_array.begin());
        gg.assign(shape_array);
    }
    Graph g;
    createCorrectGridND(shape,g);

    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        std::vector<AdjacencyType> correct;
        std::vector<AdjacencyType> testing;
        correct.insert(correct.begin(),g.adjacenciesFromVertexBegin(v),g.adjacenciesFromVertexEnd(v));
        testing.insert(testing.begin(),gg.adjacenciesFromVertexBegin(v),gg.adjacenciesFromVertexEnd(v));
        // Test endpoints
        {
            test(gg.vertexOfEdge(v,0) == g.vertexOfEdge(v,0));
            test(gg.vertexOfEdge(v,1) == g.vertexOfEdge(v,1));
        }
        {
            std::vector<size_t> vertices;
            std::vector<size_t> edges;
            vertices.insert(vertices.begin(),gg.verticesFromVertexBegin(v),gg.verticesFromVertexEnd(v));
            edges.insert(edges.begin(),gg.edgesFromVertexBegin(v),gg.edgesFromVertexEnd(v));
            for(size_t i=0; i<vertices.size(); ++i) {
                test(testing[i].vertex() == vertices[i]);
                test(testing[i].edge() == edges[i]);
            }
        }
        {
            size_t i;
            for(i=0; i<gg.numberOfEdgesFromVertex(v); ++i) {
                AdjacencyType a = gg.adjacencyFromVertex(v,i);
                test( testing[i] == a );
            }
            test(i == testing.size());
            test(std::is_sorted(testing.begin(),testing.end()));
        }
        std::sort(correct.begin(),correct.end());
        std::sort(testing.begin(),testing.end());
        test(correct == testing);
    }
}

// 2D Support functions

void createCorrectGrid(const size_t w, const size_t h, Graph& g) {
    g.assign();
    g.insertVertices(w*h);
    for(size_t x=0; x<w; ++x)
        for(size_t y=0; y<h-1; ++y)
            g.insertEdge(y+x*h,y+x*h+1); // insert vertical edge
    for(size_t x=0; x<w-1; ++x)
        for(size_t y=0; y<h; ++y)
            g.insertEdge(y+x*h,y+x*h+h); // insert vertical edge
}

void fillEdgeEndpointVector(size_t w, size_t h, std::vector<VertexPair>& vertexPairEdgeMap) {
    vertexPairEdgeMap.clear();
    for(size_t x=0; x<w; ++x)
        for(size_t y=0; y<h-1; ++y)
            vertexPairEdgeMap.push_back( {{y+x*h,y+x*h+1}}); // insert vertical edge
    for(size_t x=0; x<w-1; ++x)
        for(size_t y=0; y<h; ++y)
            vertexPairEdgeMap.push_back( {{y+x*h,y+x*h+h}}); // insert vertical edge
}

// 2D Tests

void testConstructionAndNumbers() {
    typedef andres::graph::GridGraph<2> GridGraphType;
    typedef GridGraphType::VertexCoordinate VertexCoordinate;
    typedef GridGraphType::EdgeCoordinate EdgeCoordinate;

    // Simple examples
    {
        GridGraphType gridGraph({6,5});
        VertexCoordinate originCoordinate({{0,0}});
        size_t originIndex(0);
        VertexCoordinate tenthCoordinate({{3,1}});
        size_t tenthIndex(9);
        test(originIndex == gridGraph.vertex(originCoordinate));
        test(tenthIndex == gridGraph.vertex(tenthCoordinate));
    }
    // Test initialization lists
    {
        GridGraphType gridGraph = {6,5};
        test(gridGraph.shape(0) == 6);
        test(gridGraph.shape(1) == 5);
    }
    {
        // Second edge
        GridGraphType gridGraph({6,5});
        EdgeCoordinate secondEdgeCoordinate({{1,0}}, 0);
        size_t secondEdgeIndex(1);
        test(secondEdgeIndex == gridGraph.edge(secondEdgeCoordinate));
        // edge uv joining the vertices u=(4,3) and v=(4,4) (i.e. along the second dimension)
        EdgeCoordinate uvEdgeCoordinate({{4,3}}, 1);
        size_t uvEdgeIndex(47); //
        test(uvEdgeIndex == gridGraph.edge(uvEdgeCoordinate));
    }

    {
        const GridGraphType g;
        test(!g.multipleEdgesEnabled());
        test(g.numberOfVertices() == 0);
        test(g.numberOfEdges() == 0);
    }

    for(std::size_t w = 1; w < 20; ++w) {
        for (std::size_t h=1; h< 20 ; ++h) {
            GridGraphType g(VertexCoordinate( {{h,w}}));
            const std::size_t numEdges = (h-1)*w + (w-1)*h;
            test(!g.multipleEdgesEnabled());
            test(g.numberOfVertices() == w*h);
            test(g.numberOfEdges() == numEdges);
            VertexCoordinate shape;
            for(size_t i=0; i<g.DIMENSION; ++i)
                shape[i] = g.shape(i);
            test(shape == VertexCoordinate( {{h,w}}));
            test(g.shape(0)==h && g.shape(1)==w);
            for(std::size_t x = 0; x < w; ++x) {
                for(std::size_t y = 0; y < h; ++y) {
                    std::size_t edges = 0;
                    if(x!=0) ++edges;
                    if(x!=w-1) ++edges;
                    if(y!=0) ++edges;
                    if(y!=h-1) ++edges;
                    const std::size_t index = y+x*h;
                    test(g.numberOfEdgesFromVertex(index) == edges);
                    test(g.numberOfEdgesToVertex(index) == edges);
                }
            }
        }
    }
}

// tests consistency of adjacency functions with findEdge
void testAdjacency() {
    typedef andres::graph::GridGraph<2> GridGraph;
    typedef GridGraph::AdjacencyType AdjacencyType;
    typedef typename GridGraph::VertexCoordinate VertexCoordinate;

    const size_t w = 10;
    const size_t h = 6;
    Graph g;
    const GridGraph gg(VertexCoordinate({{h,w}}));
    createCorrectGrid(w,h,g);

    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        std::vector<AdjacencyType> correct;
        std::vector<AdjacencyType> testing;
        correct.insert(correct.begin(),g.adjacenciesFromVertexBegin(v),g.adjacenciesFromVertexEnd(v));
        testing.insert(testing.begin(),gg.adjacenciesFromVertexBegin(v),gg.adjacenciesFromVertexEnd(v));
        {
            std::vector<size_t> vertices;
            std::vector<size_t> edges;
            vertices.insert(vertices.begin(),gg.verticesFromVertexBegin(v),gg.verticesFromVertexEnd(v));
            edges.insert(edges.begin(),gg.edgesFromVertexBegin(v),gg.edgesFromVertexEnd(v));
            for(size_t i=0; i<vertices.size(); ++i) {
                test(testing[i].vertex() == vertices[i]);
                test(testing[i].edge() == edges[i]);
            }
        }

        {
            size_t i;
            for(i=0; i<gg.numberOfEdgesFromVertex(v); ++i) {
                AdjacencyType a = gg.adjacencyFromVertex(v,i);
                test( testing[i] == a );
            }
            test(i == testing.size());
            test(std::is_sorted(testing.begin(),testing.end()));
        }
        std::sort(correct.begin(),correct.end());
        std::sort(testing.begin(),testing.end());
        test(correct == testing);
    }
}

template<typename G>
void testGridVertexIteratorCoordinate(const G& g) {
    typedef G Graph;
    typedef typename Graph::VertexIterator VertexIterator;
    typedef typename Graph::VertexCoordinate VertexCoordinate;
    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // verticesFromVertex
        {
            // operator[]
            for(VertexIterator it = g.verticesFromVertexBegin(v);
                it != g.verticesFromVertexEnd(v); ++it) {
                VertexCoordinate vCoord;
                const std::size_t v = *it;
                g.vertex(v,vCoord);
                VertexCoordinate itCoord;
                it.coordinate(itCoord);
                test(vCoord == itCoord);
            }
        }
    }
}


template<typename G>
void testGridVertexIterator(const G& g, const size_t pivot) {
    typedef G Graph;
    typedef typename Graph::VertexIterator VertexIterator;
    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // verticesFromVertex
        {
            // operator[]
            VertexIterator it = g.verticesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(it[j] == g.vertexFromVertex(v, j));
            }
        }
        {
            // operator*
            VertexIterator it = g.verticesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(*it == g.vertexFromVertex(v, j));
                ++it;
            }
            test(it == g.verticesFromVertexEnd(v));
        }
        // verticesToVertex
        {
            // operator[]
            VertexIterator it = g.verticesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(it[j] == g.vertexToVertex(v, j));
            }
        }
        {
            // operator*, operator++, operator==
            VertexIterator it = g.verticesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
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
        std::vector<size_t> neighbours;
        neighbours.insert(neighbours.begin(),g.verticesFromVertexBegin(pivot),g.verticesFromVertexEnd(pivot));
        {
            VertexIterator it = g.verticesFromVertexBegin(pivot);
            // operator++ (prefix)
            VertexIterator it2 = ++it;
            test(it2 == it);
            test(*it == neighbours[1]);
            // operator-- (prefix)
            it2 = --it;
            test(it2 == it);
            test(*it == neighbours[0]);
            // operator++ (postfix)
            it2 = it++;
            test(it2 != it);
            test(*it == neighbours[1]);
            test(*it2 == neighbours[0]);
            // operator-- (postfix)
            it2 = it--;
            test(it2 != it);
            test(*it == neighbours[0]);
            test(*it2 == neighbours[1]);
            // operator+=
            it2 = (it += 2);
            test(it2 == it);
            test(*it == neighbours[2]);
            // operator-=
            it2 = (it -= 2);
            test(it2 == it);
            test(*it == neighbours[0]);
        }

        // operator+, operator-
        {
            VertexIterator it = g.verticesFromVertexBegin(pivot);
            // operator+
            VertexIterator it2 = it + 2;
            test(it2 != it);
            test(*it == neighbours[0]);
            test(*it2 == neighbours[2]);
            // operator-
            VertexIterator it3 = it2 - 2;
            test(it3 != it2);
            test(*it2 == neighbours[2]);
            test(*it3 == neighbours[0]);
        }
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


template<typename G>
void testGridEdgeIterator(const G& g, const size_t pivot) {
    typedef G Graph;
    typedef typename Graph::EdgeIterator EdgeIterator;
    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // edgeFromVertex
        {
            // operator[]
            EdgeIterator it = g.edgesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(it[j] == g.edgeFromVertex(v, j));
            }
        }
        {
            // operator*
            EdgeIterator it = g.edgesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(*it == g.edgeFromVertex(v, j));
                ++it;
            }
            test(it == g.edgesFromVertexEnd(v));
        }
        // verticesToEdge
        {
            // operator[]
            EdgeIterator it = g.edgesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(it[j] == g.edgeToVertex(v, j));
            }
        }
        {
            // operator*, operator++, operator==
            EdgeIterator it = g.edgesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
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
        std::vector<size_t> neighbours;
        neighbours.insert(neighbours.begin(),g.edgesFromVertexBegin(pivot),g.edgesFromVertexEnd(pivot));
        {
            EdgeIterator it = g.edgesFromVertexBegin(pivot);
            // operator++ (prefix)
            EdgeIterator it2 = ++it;
            test(it2 == it);
            test(*it == neighbours[1]);
            // operator-- (prefix)
            it2 = --it;
            test(it2 == it);
            test(*it == neighbours[0]);
            // operator++ (postfix)
            it2 = it++;
            test(it2 != it);
            test(*it == neighbours[1]);
            test(*it2 == neighbours[0]);
            // operator-- (postfix)
            it2 = it--;
            test(it2 != it);
            test(*it == neighbours[0]);
            test(*it2 == neighbours[1]);
            // operator+=
            it2 = (it += 2);
            test(it2 == it);
            test(*it == neighbours[2]);
            // operator-=
            it2 = (it -= 2);
            test(it2 == it);
            test(*it == neighbours[0]);
        }

        // operator+, operator-
        {
            EdgeIterator it = g.edgesFromVertexBegin(pivot);
            // operator+
            EdgeIterator it2 = it + 2;
            test(it2 != it);
            test(*it == neighbours[0]);
            test(*it2 == neighbours[2]);
            // operator-
            EdgeIterator it3 = it2 - 2;
            test(it3 != it2);
            test(*it2 == neighbours[2]);
            test(*it3 == neighbours[0]);
        }
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

template<typename G>
void testGridAdjacencyIterator(const G& g, const size_t pivot) {
    typedef G GridGraph;
    typedef typename GridGraph::AdjacencyIterator AdjacencyIterator;
    typedef typename GridGraph::AdjacencyType AdjacencyType;
    // operator*, operator[]
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        // adjacencyFromVertex
        {
            // operator[]
            AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(it[j] == g.adjacencyFromVertex(v, j));
            }
        }
        {
            // operator*
            AdjacencyIterator it = g.adjacenciesFromVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(*it == g.adjacencyFromVertex(v, j));
                ++it;
            }
            test(it == g.adjacenciesFromVertexEnd(v));
        }
        // verticesToAdjacency
        {
            // operator[]
            AdjacencyIterator it = g.adjacenciesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
                test(it[j] == g.adjacencyToVertex(v, j));
            }
        }
        {
            // operator*, operator++, operator==
            AdjacencyIterator it = g.adjacenciesToVertexBegin(v);
            for(std::size_t j = 0; j < g.numberOfEdgesFromVertex(v); ++j) {
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
    {
        std::vector<AdjacencyType> neighbours;
        neighbours.insert(neighbours.begin(),g.adjacenciesFromVertexBegin(pivot),g.adjacenciesFromVertexEnd(pivot));
        {
            AdjacencyIterator it = g.adjacenciesFromVertexBegin(pivot);
            // operator++ (prefix)
            AdjacencyIterator it2 = ++it;
            test(it2 == it);
            test(*it == neighbours[1]);
            // operator-- (prefix)
            it2 = --it;
            test(it2 == it);
            test(*it == neighbours[0]);
            // operator++ (postfix)
            it2 = it++;
            test(it2 != it);
            test(*it == neighbours[1]);
            test(*it2 == neighbours[0]);
            // operator-- (postfix)
            it2 = it--;
            test(it2 != it);
            test(*it == neighbours[0]);
            test(*it2 == neighbours[1]);
            // operator+=
            it2 = (it += 2);
            test(it2 == it);
            test(*it == neighbours[2]);
            // operator-=
            it2 = (it -= 2);
            test(it2 == it);
            test(*it == neighbours[0]);
        }

        // operator+, operator-
        {
            AdjacencyIterator it = g.adjacenciesFromVertexBegin(pivot);
            // operator+
            AdjacencyIterator it2 = it + 2;
            test(it2 != it);
            test(*it == neighbours[0]);
            test(*it2 == neighbours[2]);
            // operator-
            AdjacencyIterator it3 = it2 - 2;
            test(it3 != it2);
            test(*it2 == neighbours[2]);
            test(*it3 == neighbours[0]);
        }
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

void testCopyAndAssignment() {
    // 2D
    {
        typedef andres::graph::GridGraph<> GridGraph;
        GridGraph graph0({5, 3});
        {
            GridGraph graph1(graph0); // copy
            testGraphsEquals(graph0, graph1);
        }
        {
            GridGraph graph1;
            graph1 = graph0; // assignment
            testGraphsEquals(graph0, graph1);
        }
    }

    // 3D
    {
        typedef andres::graph::GridGraph<3> GridGraph;
        GridGraph graph0({5, 3, 7});
        {
            GridGraph graph1(graph0); // copy
            testGraphsEquals(graph0, graph1);
        }
        {
            GridGraph graph1;
            graph1 = graph0; // assignment
            testGraphsEquals(graph0, graph1);
        }
    }
}

int main() {
    testConstructionAndNumbersExplicit(); //
    testConstructionAndNumbers(); // explicit test
    testAdjacency();
    testCopyAndAssignment();

    {
        typedef andres::graph::GridGraph<> GridGraph;
        typedef typename GridGraph::VertexCoordinate VertexCoordinate;
        {

            std::vector<size_t> shape = {{6,5}};
            VertexCoordinate vcPivot = {{4,3}};
            VertexCoordinate vcShape;
            std::copy(shape.begin(),shape.end(),vcShape.begin());

            testFindEdgeND<GridGraph>(shape); // comparison test

            GridGraph g(vcShape);
            const size_t pivot = g.vertex(vcPivot);
            testGridVertexIterator(g, pivot);
            testGridEdgeIterator(g, pivot);
            testGridAdjacencyIterator(g, pivot);
        }
    }

    {
        typedef andres::graph::GridGraph<4> GridGraph;
        typedef typename GridGraph::VertexCoordinate VertexCoordinate;
        {
            GridGraph g({6,4,3,5});
            const size_t pivot = g.vertex(VertexCoordinate({{1,2,2,1}}));
            testIteratorCompile(g, pivot);
        }
        {
            std::vector<size_t> shape = {{4,5,3,6}};
            testConstructionAndNumbersND<GridGraph>(shape);
            testFindEdgeND<GridGraph>(shape);
            testAdjacencyND<GridGraph>(shape);
            testOrientationsND<GridGraph>(shape);
        }
        {
            GridGraph g(VertexCoordinate({{3,2,3,2}}));
            const size_t pivot = g.vertex(VertexCoordinate({{1,1,1,1}}));
            testGridVertexIterator(g, pivot);
            testGridVertexIteratorCoordinate(g);
            testGridEdgeIterator(g, pivot);
            testGridAdjacencyIterator(g, pivot);
        }
    }
    return 0;
}
