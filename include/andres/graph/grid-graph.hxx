
#pragma once
#ifndef ANDRES_GRAPH_GRID_GRAPH_HXX
#define ANDRES_GRAPH_GRID_GRAPH_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <iterator>
#include <array>

#include "adjacency.hxx"
#include "visitor.hxx"

namespace andres {

namespace graph {

/// \class GridGraph
/// n-dimensional grid graph.
/// \par Element Indexing
/// This graph references the vertices and edges using linearly increasing
/// indices in the range
/// 0 to |V|-1 and 0 to |E|-1 respectively.
/// Specifically for the grid graph accessing of the unerlying structure
/// can be ahieved:
///  - The vertices can be specified in two ways:
///     -# Using a linearly increasing integer of type GridGraph::size_type
///        in the range 0 to \b numberOfVertices().\n
///        These indices refer to the vertices on the grid in a
///        smallest-dimension-first fashion: firstly all the vertices of
///        the first dimension are enumerated, in increasing order, then
///        along the second dimension and so on.
///     -# Using a GridGraph::VertexCoordinate, effectively an \b std::array
///        which contains the coordinates of the vertex in the grid, starting
///        with the first dimension.
///     For example, in a 2-Dimensional graph with shape 5x6, the origin
///     and the 10th vertex are:
///     \code
///        typedef andres::graph::GridGraph<size_t,2> GridGraph;
///        typedef typename GridGraph::VertexCoordinate VertexCoordinate;
///        GridGraph gridGraph({6,5});
///        VertexCoordinate originCoordinate({0,0});
///        size_t originIndex(0);
///        VertexCoordinate tenthCoordinate({3,1});
///        size_t tenthIndex(9);
///        test(originIndex == gridGraph.vertex(originCoordinate)); // OK
///        test(tenthIndex == gridGraph.vertex(tenthCoordinate)); // OK
///     \endcode
///     To convert between the two formats consult GridGraph::vertex().
///  - The edges can be specified in two ways:
///     -# Using a linearly increasing integer of type GridGraph::size_type
///        in the range 0 to \b numberOfVertices().\n
///        These indices refer to the edges on the grid in a
///        smallest-dimension-first fashion: firstly all the edges than run
///        along the first dimension (i.e.: those that
///        whose endpoints are vertices whose coordinates differ by 1 only
///        on the first dimension). Consequently enumerated are those along
///        the second dimension and so on.
///     -# By specifying the coordinates of the smallest of the two endpoints
///        and the unique dimension along which the coordinates of the
///        endpoints have a unit difference.
///        These numbers are packed in an GridGraph::EdgeVertex \c struct. In
///        this documentation the vertex with the sallest coordinates is
///        referred to as the \e pivot of the edge.
///     For example, in a 2-Dimensional graph with shape 5x6
///     \code
///        typedef andres::graph::GridGraph<size_t,2> GridGraph;
///        typedef typename GridGraph::VertexCoordinate VertexCoordinate;
///        // Second edge
///        typedef typename GridGraph::EdgeCoordinate EdgeCoordinate;
///        typedef typename GridGraph::VertexCoordinate VertexCoordinate;
///        GridGraph gridGraph({6,5});
///        EdgeCoordinate secondEdgeCoordinate({1,0},0);
///        size_t secondEdgeIndex(1);
///        test(secondEdgeIndex == gridGraph.edge(secondEdgeCoordinate)); // OK
///        // edge uv joining the vertices u=(4,3) and v=(4,4)
///        EdgeCoordinate uvEdgeCoordinate({4,3},1); size_t uvEdgeIndex(47);
///        // test(uvEdgeIndex == gridGraph.edge(uvEdgeCoordinate)); /// OK
/// \endcode To convert between the two formats consult GridGraph::edge().
/// \par Adjacency Indexing The adjacent elements for each vertex (be it
///     Vertices, Edges or Adjacencies) are guaranteed to be enumerated in
///     strictly increasing order.
///
template <
std::size_t D = 2,
    typename S = std::size_t,
    typename VISITOR = IdleGraphVisitor<S>
    >
class GridGraph {
public:
    typedef S size_type;
    typedef VISITOR Visitor;
    typedef andres::graph::Adjacency<size_type> Adjacency;

    static const size_type DIMENSION = static_cast<size_type> (D);

    typedef std::array<size_type, DIMENSION> VertexCoordinate;

    /// \class EdgeCoordinate
    /// Describes an edge as the integer index of the minimum of the
    /// two endpoints and the direction along which it is drawn.
    struct EdgeCoordinate {
        EdgeCoordinate(const VertexCoordinate&, const size_type, bool = false);
        EdgeCoordinate();
        /// The minimum of the two endpoints of the edge is specified
        /// as an integer and accessed by the \b pivot member.
        /// This specifies the \e origin of the edge, in a direction
        /// towards the positive orthant.
        /// (i.e.: in specifying the pivot of an edge, an isSmaller
        /// of false is implied.)
        VertexCoordinate pivotCoordinate_; // TODO: I noticed an inadvertently overlooked TODO: bjoern-comment intructing a name-change. In the new version it disappeared. Which name would you prefer? ( I am usually doing automatic replacements in updating the code) Recovered comment: TODO: bjoern: indent and call vertex0_
        /// The dimension along which the edge is drawn.
        size_type dimension_;
    };

    /// AdjacencyIterator
    // \cond SUPPRESS_DOXYGEN
    class AdjacencyIterator
        :   public std::iterator <
			std::random_access_iterator_tag,
			const Adjacency
        >  { // TODO: bjoern: wrap line and indent, also: does this define the typedefs of iterator_traits? // Yes. This seems to be the de-facto way to proceed. I had some doubts with overriding the typedefs in the derived classes, but I made a testcase and it seems fine!
    public:
        typedef GridGraph<DIMENSION, size_type, Visitor> GraphType;
        typedef typename AdjacencyIterator::difference_type difference_type;
        typedef typename AdjacencyIterator::pointer pointer;
        typedef typename AdjacencyIterator::reference reference;

        AdjacencyIterator();
        AdjacencyIterator(const GraphType&);
        AdjacencyIterator(const GraphType&, const size_type);
        AdjacencyIterator(const GraphType&, const size_type, const size_type);

        // increment and decrement
        AdjacencyIterator& operator+=(const difference_type);
        AdjacencyIterator& operator-=(const difference_type);
        AdjacencyIterator& operator++(); // prefix
        AdjacencyIterator& operator--(); // prefix
        AdjacencyIterator operator++(int); // postfix
        AdjacencyIterator operator--(int); // postfix
        AdjacencyIterator operator+(const difference_type) const;
        AdjacencyIterator operator-(const difference_type) const;
        difference_type operator-(const AdjacencyIterator&) const;

        // comparison
        bool operator==(const AdjacencyIterator&) const;
        bool operator!=(const AdjacencyIterator&) const;
        bool operator<(const AdjacencyIterator&) const;
        bool operator<=(const AdjacencyIterator&) const;
        bool operator>(const AdjacencyIterator&) const;
        bool operator>=(const AdjacencyIterator&) const;

        // access
        reference operator*();
        pointer operator->();
        reference operator[](const difference_type);

    protected:
        const GraphType* graph_;
        size_type vertex_;
        size_type adjacencyIndex_;
        Adjacency adjacency_;
    };

    class VertexIterator
            :       public AdjacencyIterator {
    public:
        typedef GridGraph<DIMENSION, size_type, Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef size_type value_type;
		typedef typename VertexIterator::difference_type difference_type;
        typedef value_type* pointer;
        typedef value_type& reference;

        VertexIterator();
        VertexIterator(const VertexIterator&);
        VertexIterator(const AdjacencyIterator&);
        VertexIterator(const GraphType&);
        VertexIterator(const GraphType&, const size_type);
        VertexIterator(const GraphType&, const size_type, const size_type);

        // access
        value_type operator*() const;  // TODO: we can't use reference here, since we are creating data in-situ. However, it should not be a problem? Should I change the underlying interface?
        value_type operator[](const difference_type) const;

        VertexCoordinate coordinate() const; // TODO: Ask: I forgot to note this last time. Review needed. I didn't test this
    };

    class EdgeIterator
            :       public AdjacencyIterator {
    public:
        typedef GridGraph<DIMENSION, size_type, Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef size_type value_type;
		typedef typename EdgeIterator::difference_type difference_type;

        EdgeIterator();
        EdgeIterator(const EdgeIterator&);
        EdgeIterator(const AdjacencyIterator&);
        EdgeIterator(const GraphType&);
        EdgeIterator(const GraphType&, const size_type);
        EdgeIterator(const GraphType&, const size_type, const size_type);

        // access
        value_type operator*() const;
        value_type operator[](const difference_type) const;
    };
    // \endcond

    // construction
    GridGraph(const Visitor& = Visitor());
    GridGraph(const VertexCoordinate&, const Visitor& = Visitor());
    void assign(const Visitor& = Visitor());
    void assign(const VertexCoordinate&, const Visitor& = Visitor());

    size_type shape(const size_type) const;

    // iterator access (compatible with Digraph)
    VertexIterator verticesFromVertexBegin(const size_type) const;
    VertexIterator verticesFromVertexEnd(const size_type) const;
    VertexIterator verticesToVertexBegin(const size_type) const;
    VertexIterator verticesToVertexEnd(const size_type) const;
    EdgeIterator edgesFromVertexBegin(const size_type) const;
    EdgeIterator edgesFromVertexEnd(const size_type) const;
    EdgeIterator edgesToVertexBegin(const size_type) const;
    EdgeIterator edgesToVertexEnd(const size_type) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const size_type) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const size_type) const;
    AdjacencyIterator adjacenciesToVertexBegin(const size_type) const;
    AdjacencyIterator adjacenciesToVertexEnd(const size_type) const;

    // access (compatible with Digraph)
    size_type numberOfVertices() const;
    size_type numberOfEdges() const;
    size_type numberOfEdgesFromVertex(const size_type) const;
    size_type numberOfEdgesToVertex(const size_type) const;
    size_type vertexOfEdge(const size_type, const size_type) const;
    size_type edgeFromVertex(const size_type, const size_type) const;
    size_type edgeToVertex(const size_type, const size_type) const;
    size_type vertexFromVertex(const size_type, const size_type) const;
    size_type vertexToVertex(const size_type, const size_type) const;
    Adjacency adjacencyFromVertex(const size_type, const size_type) const;
    Adjacency adjacencyToVertex(const size_type, const size_type) const;
    std::pair<bool, size_type> findEdge(const size_type, const size_type) const;
    bool multipleEdgesEnabled() const;

    size_type insertEdge(const size_type, const size_type) const;

    size_type vertex(const VertexCoordinate&) const;
    void vertex(const size_type, VertexCoordinate&) const;
    size_type edge(const EdgeCoordinate&) const;
    void edge(const size_type, EdgeCoordinate&) const;

private:
    size_type vertexFromVertex(const VertexCoordinate&, const size_type, size_type&, bool&) const;
//  void adjacencyFromVertex(const VertexCoordinate&,const size_type j,size_type& adjacentVertexIndex,size_type& adjacentEdgeIndex) const; // TODO: janis: I leave the names for clarity. Once reviewed, (I will) remove
    void adjacencyFromVertex(const VertexCoordinate&, const size_type, size_type&, size_type&) const;

    // Member variables
    VertexCoordinate shape_;
    std::array<size_type, DIMENSION> edgeIndexOffsets_;
    std::array<size_type, DIMENSION> vertexIndexOffsets_;
    std::array<VertexCoordinate, DIMENSION> edgeShapes_;
    size_type numberOfVertices_;
    const size_type& numberOfEdges_;
    Visitor visitor_;
};

/// Construct an empty grid graph.
/// \tparam S the type of the indices used. It must be an unsigned type
/// (otherwise a compilation error is caused). Defaults to \c std::size_t .
/// \tparam VISITOR a visitor class to be used. Since this class is immutable
/// appart from resizing, this is a dummy variable meant for compatibility,
/// defaulting to \c andres::graph::IdgeGraphVisitor.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::GridGraph(
    const Visitor& visitor
)
:   GridGraph( {0, 0}, visitor) // Chain-call Constructor
{}

/// Construct a grid graph with a specified shape.
///
/// \param shape the shape of the Grid Graph as an std::array.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::GridGraph(
    const VertexCoordinate& shape,
    const Visitor& visitor
)
:   numberOfEdges_ (edgeIndexOffsets_[DIMENSION - 1]) {
    assign(shape, visitor);
}

/// Clear a grid graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
///
template<std::size_t D, typename S, typename VISITOR>
inline void
GridGraph<D, S, VISITOR>::assign(
    const Visitor& visitor
) {
    VertexCoordinate shape;
    std::fill(shape.begin(), shape.end(), 0);
    assign(shape, visitor);
}

/// Clear a grid graph and assign a new shape.
///
/// \param shape the shape of the grid graph
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<std::size_t D, typename S, typename VISITOR>
inline void
GridGraph<D, S, VISITOR>::assign(
    const VertexCoordinate& shape,
    const Visitor& visitor
) {
    shape_ = shape;
    visitor_ = visitor;

    // Set vertex offsets for fast vertex indexing. // TODO: Ask about the locality of this.
    {
        size_type cumprod = 1;
        size_type i;
        for(i = 0; i < DIMENSION; ++i) {
            vertexIndexOffsets_[i] = cumprod;
            cumprod *= shape_[i];
        }
        numberOfVertices_ = cumprod;
    }
    {
        size_type edgeIndexOffset = 0;  // First edge is at offset 0
        for(size_type i = 0; i < DIMENSION; ++i) {
            VertexCoordinate& edgeShape = edgeShapes_[i];
            edgeShape = shape_; // the i-th dimension of the edges along the i-th dimension is 1 less
            if(edgeShape[i] > 0) { // If already zero, no need to reduce.
                --edgeShape[i];
            }
            {
                //
                size_type cumprod = edgeShape[0];
                for(size_type j = 1; j < DIMENSION; ++j) {
                    cumprod *= edgeShape[j];
                }
                edgeIndexOffsets_[i] = (edgeIndexOffset += cumprod);
            }
        }
    }
}

/// Get the size of a specific dimension of the grid graph.
/// \param dimension the index of the dimension index to retrieve.
/// \return the size of the specified dimension.
/// \see mapIndexToCoordinate mapVertexCoordinateToIndex
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::shape(
    const size_type dimension
) const {
    return shape_[dimension];
}

/// Get an iterator to the beginning of the sequence of vertices reachable
/// from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexEnd()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexIterator
GridGraph<D, S, VISITOR>::verticesFromVertexBegin(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices reachable from
/// a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexBegin()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexIterator
GridGraph<D, S, VISITOR>::verticesFromVertexEnd(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of vertices from which
/// a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexEnd()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexIterator
GridGraph<D, S, VISITOR>::verticesToVertexBegin(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices from which a
/// given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexBegin()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexIterator
GridGraph<D, S, VISITOR>::verticesToVertexEnd(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that originate
/// from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexEnd()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::EdgeIterator
GridGraph<D, S, VISITOR>::edgesFromVertexBegin(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that originate from
/// a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexBegin()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::EdgeIterator
GridGraph<D, S, VISITOR>::edgesFromVertexEnd(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that are
/// incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexEnd()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::EdgeIterator
GridGraph<D, S, VISITOR>::edgesToVertexBegin(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that are incident
/// to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexBegin()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::EdgeIterator
GridGraph<D, S, VISITOR>::edgesToVertexEnd(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies that
/// originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::adjacenciesFromVertexBegin(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies that originate
/// from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::adjacenciesFromVertexEnd(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies incident
/// to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::adjacenciesToVertexBegin(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies incident to
/// a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::adjacenciesToVertexEnd(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get the number of vertices.
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::numberOfVertices() const {
    return numberOfVertices_;
}

/// Get the number of edges.
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::numberOfEdges() const {
    return numberOfEdges_;
}

/// Get the number of edges that originate from a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeFromVertex()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::numberOfEdgesFromVertex(
    const size_type vertex
) const {
    assert(vertex < numberOfVertices());
    VertexCoordinate edgeCoordinate;
    this->vertex(vertex, edgeCoordinate);
    size_type numEdgesFromVertex = 0;
    for(size_type i = 0; i < DIMENSION; ++i) {
        const size_type& coordinate = edgeCoordinate[i];
        if(coordinate > 0) {
            ++numEdgesFromVertex;
        }
        if(coordinate < (shape_[i] - 1)) {
            ++numEdgesFromVertex;
        }
    }
    return numEdgesFromVertex;
}

/// Get the number of  edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::numberOfEdgesToVertex(
    const size_type vertex
) const {
    return numberOfEdgesFromVertex(vertex);
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::vertexOfEdge(
    const size_type edge,
    const size_type j
) const {
    assert(edge < numberOfEdges());
    assert(j < 2);
    EdgeCoordinate edgeCoordinate;
    this->edge(edge, edgeCoordinate);
    size_type pivotIndex = vertex(edgeCoordinate.pivotCoordinate_);
    if(j == 0) {
        return pivotIndex;
    }
    else {
        return pivotIndex + vertexIndexOffsets_[edgeCoordinate.dimension_];
    }
}

/// Get the integer index of an edge that originates from a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex)
/// - 1.
///
/// \sa numberOfEdgesFromVertex()exFi
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::edgeFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesFromVertex(vertex));
    VertexCoordinate vertexCoordinate;
    this->vertex(vertex, vertexCoordinate);
    size_type direction;
    bool isSmaller;
    vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
    if(isSmaller) {
        VertexCoordinate pivotVertexCoordinate = vertexCoordinate;
        --pivotVertexCoordinate[direction];
        return edge(EdgeCoordinate(pivotVertexCoordinate, direction));
    }
    else {
        return edge(EdgeCoordinate(vertexCoordinate, direction));
    }
}

/// Get the integer index of an edge that is incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex)
/// - 1.
///
/// \sa numberOfEdgesToVertex()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::edgeToVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return edgeFromVertex(vertex, j);
}

/// Get the integer index of a vertex from which a given vertex is reachable
/// via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and
/// numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::vertexFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    size_type direction;
    bool isSmaller;
    VertexCoordinate vertexCoordinate;
    this->vertex(vertex, vertexCoordinate);
    return vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
}


/// Get the integer index of a vertex to which a given vertex has a single
/// edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and
/// numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::vertexToVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return vertexFromVertex(vertex, j);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::Adjacency
GridGraph<D, S, VISITOR>::adjacencyFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    size_type direction;
    size_type adjacentEdgeIndex;
    size_type adjacentVertexIndex;

    {
        VertexCoordinate vertexCoordinate;
        this->vertex(vertex, vertexCoordinate);

        adjacencyFromVertex(vertexCoordinate, j, adjacentVertexIndex, adjacentEdgeIndex);
    }
    return Adjacency(adjacentVertexIndex, adjacentEdgeIndex);
}


/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::Adjacency
GridGraph<D, S, VISITOR>::adjacencyToVertex(
    const size_type vertex,
    const size_type j
) const {
    return adjacencyFromVertex(vertex, j);
}

/// Search for an edge (in constant time).
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \retval pair an \c std::pair. If an edge from \b vertex0 to \b vertex1
/// exists, \c pair.first is \c true
/// and \c pair.second is the index of such an edge. If no edge from \b
/// vertex0 to \b vertex1 exists, \c pair.first is \c false
/// and \c pair.second is undefined.
template<std::size_t D, typename S, typename VISITOR>
inline std::pair<bool, typename GridGraph<D, S, VISITOR>::size_type>
GridGraph<D, S, VISITOR>::findEdge(
    const size_type vertex0,
    const size_type vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());
    if(vertex0 < vertex1) {
        size_type offset = vertex1 - vertex0;
        size_type i;
        for(i = 0; i < DIMENSION - 1; ++i) {
            if(offset == vertexIndexOffsets_[i]) { // edge found unless offset runs across dimension
                VertexCoordinate vertexCoordinate0;
                vertex(vertex0, vertexCoordinate0);
                if(vertexCoordinate0[i] != shape_[i] - 1) { // Check for boundary case
                    const EdgeCoordinate edgeCoordinate(vertexCoordinate0, i, false);
                    const size_type edgeIndex = edge(edgeCoordinate);
                    return std::make_pair(true, edgeIndex);
                }
            }
        }
        if(offset == vertexIndexOffsets_[i]) { // edge found
            VertexCoordinate vertexCoordinate0;
            vertex(vertex0, vertexCoordinate0);
            const EdgeCoordinate edgeCoordinate(vertexCoordinate0, i, false);
            const size_type edgeIndex = edge(edgeCoordinate);
            return std::make_pair(true, edgeIndex);
        }
    }
    else {   // On expectation faster to ignore the equal case. (Assuming uniform input distribution)
        size_type offset = vertex0 - vertex1;
        size_type i;
        for(i = 0; i < DIMENSION - 1; ++i) {
            if(offset == vertexIndexOffsets_[i]) { // edge found unless offset runs across dimensions
                VertexCoordinate vertexCoordinate1;
                vertex(vertex1, vertexCoordinate1);
                if(vertexCoordinate1[i] != shape_[i] - 1) { // Check for boundary case
                    const EdgeCoordinate edgeCoordinate(vertexCoordinate1, i, false);
                    const size_type edgeIndex = edge(edgeCoordinate);
                    return std::make_pair(true, edgeIndex);
                }
            }
        }
        if(offset == vertexIndexOffsets_[i]) { // edge found unless offset runs across dimensions
            VertexCoordinate vertexCoordinate1;
            vertex(vertex1, vertexCoordinate1);
            const EdgeCoordinate edgeCoordinate(vertexCoordinate1, i, false);
            const size_type edgeIndex = edge(edgeCoordinate);
            return std::make_pair(true, edgeIndex);
        }
    }
    return std::make_pair(false, 0);
}

/// Indicate if multiple edges are enabled.
///
/// \return false
///
template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::multipleEdgesEnabled() const {
    return false;
}

/// Returns the edge between the specified vertices.
///
/// \param vertexIndex0 Integer index of the first vertex in the edge.
/// \param vertexIndex1 Integer index of the second vertex in the edge.
/// \return Integer index of the newly inserted edge.
/// \throw runtime_error If the edge does not exist.
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::insertEdge(
    const size_type vertexIndex0,
    const size_type vertexIndex1
) const {
    assert(vertexIndex0 < numberOfVertices()); 
    assert(vertexIndex1 < numberOfVertices()); 
    
    std::pair<bool, std::size_t> p = findEdge(vertexIndex0, vertexIndex1);
    if(p.first == true) {
        return p.second;
    } else {
        throw std::runtime_error("Edge not found.");
    }
}

/// Retrieve the specified vertex of the graph.
/// \param[in] vertexIndex the integer index of the requested vertex
/// \param[out] vertexCoordinate The coordinates of the vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<std::size_t D, typename S, typename VISITOR>
void
GridGraph<D, S, VISITOR>::vertex(
    size_type vertexIndex,
    VertexCoordinate& vertexCoordinate
) const {
    assert(vertexIndex < numberOfVertices_);
    size_type i;
    for(i = 0; i < DIMENSION - 1; ++i) {
        vertexCoordinate[i] = vertexIndex % shape_[i];
        vertexIndex = vertexIndex / shape_[i];
    }
    vertexCoordinate[i] = vertexIndex;
}

/// Retrieve the specified vertex of the graph.
/// \param vertexCoordinate the integer index of the requested vertex
/// \return The integer index of the specified vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<std::size_t D, typename S, typename VISITOR>
typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::vertex(
    const VertexCoordinate& vertexCoordinate
) const {
    size_type index = vertexCoordinate[DIMENSION - 1];
    for(size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * shape_[i - 1] + vertexCoordinate[i - 1];
    }
    return index;
}

/// Retrieve the specified edge of the graph.
/// \param edgeCoordinate the coordinates of the minimum endpoint (\e pivot)
/// and the direction of the requested edge.
/// \return The integer index of the specified edge.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa hasEdge()
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::edge(
    const EdgeCoordinate& edgeCoordinate
) const {
    assert(edgeCoordinate.dimension_ < DIMENSION);
    const size_type& dimension = edgeCoordinate.dimension_;
    const VertexCoordinate& pivotCoordinate = edgeCoordinate.pivotCoordinate_;
    const VertexCoordinate& edgeShape = edgeShapes_[dimension];

    // size_type index = mapCoordinateToIndex(edgeCoordinate,edgeShape);
    size_type index = pivotCoordinate[DIMENSION - 1];
    for(size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * edgeShape[i - 1] + pivotCoordinate[i - 1];
    }
    if(dimension > 0) {
        const size_type& offset = edgeIndexOffsets_[dimension - 1];
        index += offset;
    }
    return index;
}

/// Retrieve the specified edge of the graph.
/// \param[in] edgeIndex the integer index of the requested edge.
/// \param[out] edgeCoordinate a GridGraph::EdgeCoordinate instance. \c
/// edgeCoordinate.pivot is the integer index
///  of the minimum of the two edge endpoints; \p edgeCoordinate.direction
///  is the dimension along which
///  the edge is drawn, with an assumed positive isSmaller.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa GridGraph::bool
template<std::size_t D, typename S, typename VISITOR>
inline void
GridGraph<D, S, VISITOR>::edge(
    size_type edgeIndex,
    EdgeCoordinate& edgeCoordinate
) const {
    // WARNING: If the assertion does not hold, the code will scan until an unlikely condition.
    assert(edgeIndex < numberOfEdges());
    size_type& direction = edgeCoordinate.dimension_;
    // Find the direction as the last edge offset:
    for(direction = 0; edgeIndex >= edgeIndexOffsets_[direction]; ++direction);
    if(direction > 0) { // Not needed, but saves one memory lookup.
        const size_type& offset = edgeIndexOffsets_[direction - 1];
        edgeIndex -= offset;
    }
    const VertexCoordinate& edgeShape = edgeShapes_[direction];
    // mapEdgeIndexToCoordinate
    {
        VertexCoordinate& pivotCoordinate = edgeCoordinate.pivotCoordinate_;
        size_type i;
        for(i = 0; i < DIMENSION - 1; ++i) {
            pivotCoordinate[i] = edgeIndex % edgeShape[i];
            edgeIndex = edgeIndex / edgeShape[i];
        }
        pivotCoordinate[i] = edgeIndex;
    }
}


/// \internal
/// \brief Retrieve the direction and isSmaller of the <b>j</b>-th edge of
/// a vertex
/// \param vertexCoordinate the integer index of the vertex from which the
/// edges are originating
/// \param j the index within the set of adacent edges of the specified vertex
/// \param[out] direction the direction of the specified edge
/// \param[out] isSmaller the isSmaller of the specified edge
/// \retval found If the edge is found, a value of \c true is returned.
/// \remark This function attempts a fast imlementation that tries not
/// to waste any comparison operations and is meant for use as a building
/// block of the class.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa directionOfEdgeFromVertex, edgeFromVertex
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::size_type
GridGraph<D, S, VISITOR>::vertexFromVertex(
    const VertexCoordinate& vertexCoordinate,
    const size_type j,
    size_type& direction,
    bool& isSmaller
) const {
    assert(vertex(vertexCoordinate) < numberOfVertices());
    VertexCoordinate& modifieableVertexCoordinate = const_cast<VertexCoordinate&>(vertexCoordinate); // TODO: Ask: Is this ok to avoid copies?
    size_type cur_j = 0;
    for(size_type i = DIMENSION; i > 0; --i) {
        if(vertexCoordinate[i - 1] > 0) { // edge exists
            if(cur_j == j) {
                direction = i - 1;
                isSmaller = true;
                --modifieableVertexCoordinate[i - 1]; // Increase, avoiding copies
                const size_type vertexIndex = vertex(modifieableVertexCoordinate);
                ++modifieableVertexCoordinate[i - 1]; // Restore
                return vertexIndex;
            }
            ++cur_j;
        }
    }
    for(size_type i = 0; i < DIMENSION; ++i) {
        if(vertexCoordinate[i] < shape_[i] - 1) { // edge exists
            if(cur_j == j) {
                direction = i;
                isSmaller = false;
                ++modifieableVertexCoordinate[i]; // Increase, avoiding copies
                const size_type vertexIndex = vertex(modifieableVertexCoordinate);
                --modifieableVertexCoordinate[i]; // Restore
                return vertexIndex;
            }
            ++cur_j;
        }
    }
    throw std::out_of_range("vertex neighbor index out of range.");
}

template<std::size_t D, typename S, typename VISITOR>
inline void
GridGraph<D, S, VISITOR>::adjacencyFromVertex(
    const VertexCoordinate& vertexCoordinate,
    const size_type j,
    size_type& adjacentVertexIndex,
    size_type& adjacentEdgeIndex
) const {
    assert(vertex(vertexCoordinate) < numberOfVertices());
    size_type direction;
    bool isSmaller;
    adjacentVertexIndex = vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
    if(isSmaller) {
        VertexCoordinate pivotVertexCoordinate = vertexCoordinate;
        --pivotVertexCoordinate[direction];
        adjacentEdgeIndex = edge(EdgeCoordinate(pivotVertexCoordinate, direction));
    }
    else {
        adjacentEdgeIndex = edge(EdgeCoordinate(vertexCoordinate, direction));
    }

}

/// Initialize an edge coordinate.
/// \param pivotCoordinate coordinate of the reference vertex.
/// \param dimension dimension along which the edge is drawn.
/// \param isSmaller relative position of specified endpoint.\n
/// A value of \c true results in an object corresponding to the edge
/// of which the smallest endpoint has a GridGraph::VertexCoordinate
/// equal to \b pivotCoordinate.\n
/// Similarly, a value of \c false results in an object corresponding
/// to the edge of which the largest endpoint has a
/// GridGraph::VertexCoordinate equal to \b pivotCoordinate.
template<std::size_t D, typename S, typename VISITOR>
inline GridGraph<D, S, VISITOR>::EdgeCoordinate::EdgeCoordinate(
    const VertexCoordinate& pivotCoordinate,
    const size_type dimension,
    const bool isSmaller
)
:   pivotCoordinate_(pivotCoordinate),
        dimension_(dimension) {
    if(isSmaller) {
        assert(pivotCoordinate_[dimension] > 0);
        --pivotCoordinate_[dimension];
    }
}

/// Default non-initializing constructor
/// \warning this constructor will \b NOT set the pivotCoordinate_ membr
/// variable to zero.
template<std::size_t D, typename S, typename VISITOR>
inline GridGraph<D, S, VISITOR>::EdgeCoordinate::EdgeCoordinate() {
}


// \cond SUPPRESS_DOXYGEN

// implementation of AdjacencyIterator

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::AdjacencyIterator::AdjacencyIterator()
:   vertex_(0),
        adjacencyIndex_(0),
        adjacency_()
{}
template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph
)
:   graph_(&graph),
        vertex_(0),
        adjacencyIndex_(0),
        adjacency_()
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const size_type vertex
)
:   graph_(&graph),
        vertex_(vertex),
        adjacencyIndex_(0),
        adjacency_() {
    assert(vertex < graph.numberOfVertices());
}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
:   graph_(&graph),
        vertex_(vertex),
        adjacencyIndex_(adjacencyIndex),
        adjacency_() {
    assert(vertex < graph.numberOfVertices());
    assert(adjacencyIndex <= graph.numberOfVertices());
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator&
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator+=(
    const difference_type d
) {
    adjacencyIndex_ += d;
    return *this;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator&
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator-=(
    const difference_type d
) {
    adjacencyIndex_ -= d;
    return *this;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator&
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator++() {
    ++adjacencyIndex_;
    return *this;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator&
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator--() {
    --adjacencyIndex_;
    return *this;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator++(int) {
    AdjacencyIterator copy = *this;
    ++adjacencyIndex_;
    return copy;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator--(int) {
    AdjacencyIterator copy = *this;
    --adjacencyIndex_;
    return copy;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator+(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy += d;
    return copy;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator-(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy -= d;
    return copy;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator::difference_type
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator-(
    const AdjacencyIterator& adjacencyIterator
) const {
    return adjacencyIndex_ - adjacencyIterator.adjacencyIndex_;
}

template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator==(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ == other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator!=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ != other.adjacencyIndex_
           || vertex_ != other.vertex_
           || graph_ != other.graph_;
}

template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator<(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ < other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator<=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ <= other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator>(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ > other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<std::size_t D, typename S, typename VISITOR>
inline bool
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator>=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ >= other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

// Since the GridGraph has no backing storage for each of the Adjacency
// objects handled in the class, the AdjacencyIterator is limited to only
// one adjacency object per iterator instance.
// This should be sufficient for most normal usage scenarios, and is supported
// by the very lightweight copying of the Adjacency object itself.
// Note, however, that if you store a reference to the pointed object,
// then advance the iterator and subsequently dereference it again,
// the new reference will refer to the very same same adjacency object,
// unique for the AdjacencyIterator instance.
// This operation will therefore silently update all previous references to
// the adjacency object of the iterator.
template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator::reference
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator*() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return adjacency_;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator::pointer
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator->() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return &adjacency_;
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::AdjacencyIterator::reference
GridGraph<D, S, VISITOR>::AdjacencyIterator::operator[](
    const difference_type j
) {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_ + j);
    return adjacency_;
}

// implementation of VertexIterator

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::VertexIterator::VertexIterator()
:   Base()
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph
)
:   Base(graph)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const size_type vertex
)
:   Base(graph, vertex)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
:   Base(graph, vertex, adjacencyIndex)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::VertexIterator::VertexIterator(
    const VertexIterator& it
)
:   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::VertexIterator::VertexIterator(
    const AdjacencyIterator& it
)
:   Base(it)
{}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexIterator::value_type
GridGraph<D, S, VISITOR>::VertexIterator::operator*() const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexIterator::value_type
GridGraph<D, S, VISITOR>::VertexIterator::operator[](
    const difference_type j
) const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::VertexCoordinate
GridGraph<D, S, VISITOR>::VertexIterator::coordinate() const {
    const size_type opposite = Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_);
    return Base::graph_->vertex(opposite);
}

// implementation of EdgeIterator

template<std::size_t D, typename S, typename VISITOR> // TODO: Ask: check how this can instantiate the Graph. Shouldn't the default constructor be deleted?
inline
GridGraph<D, S, VISITOR>::EdgeIterator::EdgeIterator()
:   Base()
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph
)
:   Base(graph)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const size_type vertex
)
:   Base(graph, vertex)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
:   Base(graph, vertex, adjacencyIndex)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::EdgeIterator::EdgeIterator(
    const EdgeIterator& it
)
:   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<std::size_t D, typename S, typename VISITOR>
inline
GridGraph<D, S, VISITOR>::EdgeIterator::EdgeIterator(
    const AdjacencyIterator& it
)
:   Base(it)
{}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::EdgeIterator::value_type
GridGraph<D, S, VISITOR>::EdgeIterator::operator*() const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<std::size_t D, typename S, typename VISITOR>
inline typename GridGraph<D, S, VISITOR>::EdgeIterator::value_type
GridGraph<D, S, VISITOR>::EdgeIterator::operator[](
    const difference_type j
) const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

// \endcond


} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_GRID_GRAPH_HXX


