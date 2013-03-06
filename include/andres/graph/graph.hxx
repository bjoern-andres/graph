/// \mainpage
/// andres::graph -- Graphs with Simple Indexing of Vertices and Edges in C++ 
///
/// Copyright (c) 2013 by Bjoern Andres, http://www.andres.sc
///
/// \section section_abstract Short Description
/// This C++ header file implements directed and undirected graphs as adjacency
/// lists. Unlike in other implementations such as the boost graph library, 
/// vertices and edges are always indexed by contiguous integers. This indexing
/// simplifies the implementation of algorithms for static graphs. In dynamic 
/// settings where vertices and edges are removed from a graph, indices of 
/// vertices and edges can change. These changes can be followed, if necessary, 
/// by means of a simple visitor.
///
/// \section section_features Features
/// - Directed and undirected graphs, implemented as adjacency lists
/// - Access to vertices and edges by contiguous integer indices
/// - Access to neighboring vertices and incident edges by STL-compliant random access iterators
/// - Insertion and removal of vertices and edges 
/// - Mutiple edges, which are disabled by default, can be enabled
/// - Visitors that follow changes of vertex and edge indices 
/// - Algorithms
///   - Minimum multicuts by interger programming, using Cplex or Gurobi 
///   - Connected components by breadth-first search and disjoint sets 
///   - Shortest paths (SSSP, SPSP) in weighted and unweighted graphs.
///   . 
/// 
/// \section section_license License
///
/// Copyright (c) 2013 by Bjoern Andres.
/// 
/// This software was developed by Bjoern Andres.
/// Enquiries shall be directed to bjoern@andres.sc.
///
/// All advertising materials mentioning features or use of this software must
/// display the following acknowledgement: ``This product includes andres::graph
/// developed by Bjoern Andres. Please direct enquiries concerning andres::graph
/// to bjoern@andres.sc''.
///
/// Redistribution and use in source and binary forms, with or without 
/// modification, are permitted provided that the following conditions are met:
///
/// - Redistributions of source code must retain the above copyright notice,
///   this list of conditions and the following disclaimer.
/// - Redistributions in binary form must reproduce the above copyright notice, 
///   this list of conditions and the following disclaimer in the documentation
///   and/or other materials provided with the distribution.
/// - All advertising materials mentioning features or use of this software must 
///   display the following acknowledgement: ``This product includes 
///   andres::graph developed by Bjoern Andres. Please direct enquiries 
///   concerning andres::graph to bjoern@andres.sc''.
/// - The name of the author must not be used to endorse or promote products 
///   derived from this software without specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
/// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
/// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
/// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
/// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
/// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
/// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
/// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
/// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/// 
/// \section section_todo Suggested bug fixes, changes and improvements
/// - Explicit unit test of visitors
/// - Template parameters to enable/disable multiple edges and digons
/// - Tutorial
///
#pragma once
#ifndef ANDRES_GRAPH_HXX
#define ANDRES_GRAPH_HXX

#include <cassert>
#include <iterator> // std::random_access_iterator
#include <vector>
#include <set> 
#include <iostream>
#include <utility> // std::pair

#include "andres/random-access-set.hxx"

/// The public API.
namespace andres {

/// Graphs and Algorithms that Operate on Graphs.
namespace graph {

/// The adjacency of a vertex consists of a vertex and a connecting edge.
template<class T = size_t>
class Adjacency {
public:
    typedef T Value;

    Adjacency(const Value, const Value);
    Value vertex() const;
    Value& vertex();
    Value edge() const;
    Value& edge();
    bool operator<(const Adjacency<Value>&) const;
    bool operator<=(const Adjacency<Value>&) const;
    bool operator>(const Adjacency<Value>&) const;
    bool operator>=(const Adjacency<Value>&) const;
    bool operator==(const Adjacency<Value>&) const;
    bool operator!=(const Adjacency<Value>&) const;

private:
    Value vertex_; 
    Value edge_;
};

/// An entire graph.
template<class T = size_t>
struct DefaultSubgraphMask {
    typedef T Value;

    bool vertex(const Value v) const { return true; }
    bool edge(const Value e) const { return true; }
};

/// Return 1 for every edge.
template<class T>
struct UnitEdgeWeightIterator {
    typedef ptrdiff_t difference_type;
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::random_access_iterator_tag iterator_category;

    value_type operator[](const size_t) const
        { return 1; }
};

// \cond SUPPRESS_DOXYGEN
namespace graph_detail {

template<bool DIRECTED, class T = size_t>
class Edge {
public:
    typedef T Value;

    Edge(const Value, const Value);
    Value operator[](const size_t) const;
    Value& operator[](const size_t);

private:
    Value vertexIndices_[2];
};

typedef RandomAccessSet<Adjacency<> > Adjacencies;

template<bool T>
class IteratorHelper
:   public Adjacencies::const_iterator 
{
private:
    typedef typename Adjacencies::const_iterator Base;

public:
    typedef typename Base::iterator_category iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::difference_type difference_type;

    // construction and assignment
    IteratorHelper();
    IteratorHelper(const Base&);
    IteratorHelper(const IteratorHelper<T>&);
    IteratorHelper operator=(const Base&);
    IteratorHelper operator=(const IteratorHelper<T>&);

    // increment and decrement
    IteratorHelper<T>& operator+=(const difference_type);
    IteratorHelper<T>& operator-=(const difference_type);
    IteratorHelper<T>& operator++(); // prefix
    IteratorHelper<T>& operator--(); // prefix
    IteratorHelper<T> operator++(int); // postfix
    IteratorHelper<T> operator--(int); // postfix
    IteratorHelper<T> operator+(const difference_type) const;
    IteratorHelper<T> operator-(const difference_type) const;
    difference_type operator-(const IteratorHelper<T>&) const;

    // access
    size_t operator*() const;
    size_t operator[](const size_t j) const;
};

typedef IteratorHelper<true> VertexIterator;
typedef IteratorHelper<false> EdgeIterator;

} // namespace graph_detail
// \endcond

/// Visitors can be used to follow the indices of vertices and edges.
///
/// These indices change due to the insetion and removal of vertices and edges.
///
struct IdleGraphVisitor {
    IdleGraphVisitor() {}
    void insertVertex(const size_t a) const {}
    void insertVertices(const size_t a, const size_t n) const {}
    void eraseVertex(const size_t a) const {}
    void relabelVertex(const size_t a, const size_t b) const {}
    void insertEdge(const size_t a) const {}
    void eraseEdge(const size_t a) const {}
    void relabelEdge(const size_t a, const size_t b) const {}
};

/// Visitors can be used to follow the indices of vertices and edges.
///
/// These indices change due to the insetion and removal of vertices and edges.
///
struct VerboseGraphVisitor {
    VerboseGraphVisitor() {}
    void insertVertex(const size_t a) const 
        { std::cout << "inserting vertex " << a << std::endl; }
    void insertVertices(const size_t a, const size_t n) const 
        { std::cout << "inserting " << n << " vertices, starting from index " << a << std::endl; }
    void eraseVertex(const size_t a) const 
        { std::cout << "removing vertex " << a << std::endl; }
    void relabelVertex(const size_t a, const size_t b) const 
        { std::cout << "relabeling vertex " << a << ". new label is " << b << std::endl; }
    void insertEdge(const size_t a) const 
        { std::cout << "inserting edge " << a << std::endl; }
    void eraseEdge(const size_t a) const 
        { std::cout << "removing edge " << a << std::endl; }
    void relabelEdge(const size_t a, const size_t b) const 
        { std::cout << "relabeling edge " << a << ". new label is " << b << std::endl; }
};
   
/// Undirected graph, implemented as an adjacency list.
template<typename VISITOR = IdleGraphVisitor>
class Graph {
public: 
    typedef VISITOR Visitor;
    typedef graph_detail::VertexIterator VertexIterator;
    typedef graph_detail::EdgeIterator EdgeIterator;
    typedef graph_detail::Adjacencies::const_iterator AdjacencyIterator;

    // construction
    Graph(const Visitor& = Visitor());
    Graph(const size_t, const Visitor& = Visitor());
    void reserveVertices(const size_t);
    void reserveEdges(const size_t);

    // iterator access (compatible with Digraph)
    VertexIterator verticesFromVertexBegin(const size_t) const;
    VertexIterator verticesFromVertexEnd(const size_t) const;
    VertexIterator verticesToVertexBegin(const size_t) const;
    VertexIterator verticesToVertexEnd(const size_t) const;
    EdgeIterator edgesFromVertexBegin(const size_t) const;
    EdgeIterator edgesFromVertexEnd(const size_t) const;
    EdgeIterator edgesToVertexBegin(const size_t) const;
    EdgeIterator edgesToVertexEnd(const size_t) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const size_t) const; 
    AdjacencyIterator adjacenciesFromVertexEnd(const size_t) const;
    AdjacencyIterator adjacenciesToVertexBegin(const size_t) const;
    AdjacencyIterator adjacenciesToVertexEnd(const size_t) const;

    // access (compatible with Digraph)
    size_t numberOfVertices() const;
    size_t numberOfEdges() const;
    size_t numberOfEdgesFromVertex(const size_t) const; 
    size_t numberOfEdgesToVertex(const size_t) const; 
    size_t vertexOfEdge(const size_t, const size_t) const;
    size_t edgeFromVertex(const size_t, const size_t) const;
    size_t edgeToVertex(const size_t, const size_t) const; 
    size_t vertexFromVertex(const size_t, const size_t) const; 
    size_t vertexToVertex(const size_t, const size_t) const; 
    const Adjacency<>& adjacencyFromVertex(const size_t, const size_t) const;  
    const Adjacency<>& adjacencyToVertex(const size_t, const size_t) const;  
    std::pair<bool, size_t> findEdge(const size_t, const size_t) const;
    bool multipleEdgesEnabled() const;

    // manipulation
    size_t insertVertex();
    size_t insertVertices(const size_t);
    size_t insertEdge(const size_t, const size_t);
    void eraseVertex(const size_t);
    void eraseEdge(const size_t);
    bool& multipleEdgesEnabled();

private:
    typedef Adjacency<> Adjacency;
    typedef graph_detail::Adjacencies Vertex;
    typedef graph_detail::Edge<false> Edge;

    void insertAdjacenciesForEdge(const size_t);
    void eraseAdjacenciesForEdge(const size_t);

    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    bool multipleEdgesEnabled_;
    Visitor visitor_;
};

/// Directed graph, implemented as an adjacency list.
template<typename VISITOR = IdleGraphVisitor>
class Digraph {
public: 
    typedef VISITOR Visitor;
    typedef graph_detail::VertexIterator VertexIterator;
    typedef graph_detail::EdgeIterator EdgeIterator;
    typedef graph_detail::Adjacencies::const_iterator AdjacencyIterator;

    // construction
    Digraph(const Visitor& = Visitor());
    Digraph(const size_t, const Visitor& = Visitor());
    void reserveVertices(const size_t);
    void reserveEdges(const size_t);

    // iterator access 
    VertexIterator verticesFromVertexBegin(const size_t) const;
    VertexIterator verticesFromVertexEnd(const size_t) const;
    VertexIterator verticesToVertexBegin(const size_t) const;
    VertexIterator verticesToVertexEnd(const size_t) const;
    EdgeIterator edgesFromVertexBegin(const size_t) const;
    EdgeIterator edgesFromVertexEnd(const size_t) const;
    EdgeIterator edgesToVertexBegin(const size_t) const;
    EdgeIterator edgesToVertexEnd(const size_t) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const size_t) const; 
    AdjacencyIterator adjacenciesFromVertexEnd(const size_t) const;
    AdjacencyIterator adjacenciesToVertexBegin(const size_t) const;
    AdjacencyIterator adjacenciesToVertexEnd(const size_t) const;

    // access 
    size_t numberOfVertices() const;
    size_t numberOfEdges() const;
    size_t numberOfEdgesFromVertex(const size_t) const; 
    size_t numberOfEdgesToVertex(const size_t) const; 
    size_t vertexOfEdge(const size_t, const size_t) const;
    size_t edgeFromVertex(const size_t, const size_t) const;
    size_t edgeToVertex(const size_t, const size_t) const; 
    size_t vertexFromVertex(const size_t, const size_t) const; 
    size_t vertexToVertex(const size_t, const size_t) const; 
    const Adjacency<>& adjacencyFromVertex(const size_t, const size_t) const;  
    const Adjacency<>& adjacencyToVertex(const size_t, const size_t) const;  
    std::pair<bool, size_t> findEdge(const size_t, const size_t) const;
    bool multipleEdgesEnabled() const;

    // manipulation
    size_t insertVertex();
    size_t insertVertices(const size_t);
    size_t insertEdge(const size_t, const size_t);
    void eraseVertex(const size_t);
    void eraseEdge(const size_t);
    bool& multipleEdgesEnabled();

private:
    typedef Adjacency<> Adjacency;
    typedef graph_detail::Adjacencies Adjacencies;
    struct Vertex {
        Vertex() 
            : from_(), to_() 
            {}
        Adjacencies from_;
        Adjacencies to_;
    };
    typedef graph_detail::Edge<true> Edge;

    void insertAdjacenciesForEdge(const size_t);
    void eraseAdjacenciesForEdge(const size_t);

    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    bool multipleEdgesEnabled_;
    Visitor visitor_;
};

// \cond SUPPRESS_DOXYGEN
namespace graph_detail {

// implementation of Edge

template<bool DIRECTED, class T>
inline 
Edge<DIRECTED, T>::Edge(
    const Value v0, 
    const Value v1
) {
    if(DIRECTED) { // evaluated at compile time
        vertexIndices_[0] = v0; 
        vertexIndices_[1] = v1;
    }
    else {
        if(v0 <= v1) {
            vertexIndices_[0] = v0; 
            vertexIndices_[1] = v1;
        } 
        else {
            vertexIndices_[0] = v1; 
            vertexIndices_[1] = v0;
        }
    }
}

template<bool DIRECTED, class T>
inline typename Edge<DIRECTED, T>::Value
Edge<DIRECTED, T>::operator[]
(
    const size_t j
) const {
    assert(j < 2);

    return vertexIndices_[j];
}

template<bool DIRECTED, class T>
inline typename Edge<DIRECTED, T>::Value&
Edge<DIRECTED, T>::operator[](
    const size_t j
) {
    assert(j < 2);

    return vertexIndices_[j];
}

// implementation of IteratorHelper

template<bool T>
inline 
IteratorHelper<T>::IteratorHelper()
:   Base()
{}

template<bool T>
inline 
IteratorHelper<T>::IteratorHelper(
    const Base& it
)
:   Base(it)
{}

template<bool T>
inline 
IteratorHelper<T>::IteratorHelper(
    const IteratorHelper<T>& it
)
:   Base(it)
{}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator=(
    const Base& it
) { 
    Base::operator=(it); 
    return *this; 
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator=(
    const IteratorHelper<T>& it
) { 
    Base::operator=(it); 
    return *this; 
}

template<bool T>
inline size_t 
IteratorHelper<T>::operator*() const { 
    if(T) { // evaluated at compile time
        return Base::operator*().vertex(); 
    }
    else {
        return Base::operator*().edge(); 
    }
}

template<bool T>
inline size_t 
IteratorHelper<T>::operator[](
    const size_t j
) const {
    if(T) { // evaluated at compile time
        return Base::operator[](j).vertex(); 
    }    
    else {
        return Base::operator[](j).edge(); 
    }
}

template<bool T>
inline IteratorHelper<T>& 
IteratorHelper<T>::operator+=(
    const difference_type d
) {
    Base::operator+=(d);
    return *this;
}

template<bool T>
inline IteratorHelper<T>& 
IteratorHelper<T>::operator-=(
    const difference_type d
) {
    Base::operator-=(d);
    return *this;
}

template<bool T>
inline IteratorHelper<T>& 
IteratorHelper<T>::operator++() { // prefix
    Base::operator++();
    return *this;
}

template<bool T>
inline IteratorHelper<T>& 
IteratorHelper<T>::operator--() { // prefix
    Base::operator--();
    return *this;
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator++(int) { // postfix
    return Base::operator++(int());
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator--(int) { // postfix
    return Base::operator--(int());
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator+(
    const difference_type d
) const {
    return Base::operator+(d);
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator-(
    const difference_type d
) const {
    return Base::operator-(d);
}

template<bool T>
inline typename IteratorHelper<T>::difference_type
IteratorHelper<T>::operator-(
    const IteratorHelper<T>& other
) const {
    return Base::operator-(other);
}

} // namespace graph_detail
// \endcond

// implementation of Adjacency

/// Construct an adjacency.
///
/// \param vertex Vertex.
/// \param edge Edge.
///
template<class T>
inline 
Adjacency<T>::Adjacency(
    const Value vertex, 
    const Value edge
) 
:   vertex_(vertex), 
    edge_(edge)
{}

/// Access the vertex.
///
template<class T>
inline typename Adjacency<T>::Value  
Adjacency<T>::vertex() const { 
    return vertex_; 
}

/// Access the vertex.
///
template<class T>
inline typename Adjacency<T>::Value& 
Adjacency<T>::vertex() {
    return vertex_; 
}

/// Access the edge.
///
template<class T>
inline typename Adjacency<T>::Value
Adjacency<T>::edge() const { 
    return edge_; 
}

/// Access the edge.
///
template<class T>
inline typename Adjacency<T>::Value& 
Adjacency<T>::edge() { 
    return edge_; 
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool 
Adjacency<T>::operator<(
    const Adjacency<T>& in
) const { 
    if(vertex_ < in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ < in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool 
Adjacency<T>::operator<=(
    const Adjacency<T>& in
) const {
    if(vertex_ < in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ <= in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool 
Adjacency<T>::operator>(
    const Adjacency<T>& in
) const { 
    if(vertex_ > in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ > in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool 
Adjacency<T>::operator>=(
    const Adjacency<T>& in
) const { 
    if(vertex_ > in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ >= in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are equal if both the vertex and the edge are equal.
///
template<class T>
inline bool 
Adjacency<T>::operator==(
    const Adjacency<T>& in
) const { 
    return vertex_ == in.vertex_ && edge_ == in.edge_; 
}

/// Adjacencies are unequal if either the vertex or the edge are unqual.
///
template<class T>
inline bool 
Adjacency<T>::operator!=(
    const Adjacency<T>& in
) const { 
    return !(*this == in); 
}

// implementation of Graph

/// Construct an undirected graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline 
Graph<VISITOR>::Graph(
    const Visitor& visitor
)
:   vertices_(),
    edges_(),
    multipleEdgesEnabled_(false),
    visitor_(visitor)
{}

/// Construct an undirected graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline 
Graph<VISITOR>::Graph(
    const size_t numberOfVertices,
    const Visitor& visitor
)
:   vertices_(numberOfVertices),
    edges_(),
    multipleEdgesEnabled_(false),
    visitor_(visitor)
{
    visitor_.insertVertices(0, numberOfVertices);
}

/// Get the number of vertices.
///
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::numberOfVertices() const { 
    return vertices_.size(); 
}

/// Get the number of edges.
///
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::numberOfEdges() const { 
    return edges_.size(); 
}

/// Get the number of edges that originate from a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeFromVertex()
///
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::numberOfEdgesFromVertex(
    const size_t vertex
) const { 
    return vertices_[vertex].size();
}

/// Get the number of edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::numberOfEdgesToVertex(
    const size_t vertex
) const { 
    return vertices_[vertex].size();
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<typename VISITOR>
inline size_t
Graph<VISITOR>::vertexOfEdge(
    const size_t edge,
    const size_t j
) const {
    assert(j < 2);

    return edges_[edge][j];
}

/// Get the integer index of an edge that originates from a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline size_t
Graph<VISITOR>::edgeFromVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex][j].edge();
}

/// Get the integer index of an edge that is incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesToVertex()
///
template<typename VISITOR>
inline size_t
Graph<VISITOR>::edgeToVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex][j].edge();
}

/// Get the integer index of a vertex reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex() 
///
template<typename VISITOR>
inline size_t
Graph<VISITOR>::vertexFromVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex][j].vertex();
}

/// Get the integer index of a vertex from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex() 
///
template<typename VISITOR>
inline size_t
Graph<VISITOR>::vertexToVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex][j].vertex();
}

/// Insert an additional vertex.
///
/// \return Integer index of the newly inserted vertex.
///
/// \sa insertVertices()
///
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::insertVertex() {
    vertices_.push_back(Vertex());
    visitor_.insertVertex(vertices_.size() - 1);
    return vertices_.size() - 1;
}

/// Insert additional vertices.
///
/// \param number Number of new vertices to be inserted.
/// \return Integer index of the first newly inserted vertex.
///
/// \sa insertVertex()
///
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::insertVertices(
    const size_t number
) {
    size_t position = vertices_.size();
    vertices_.insert(vertices_.end(), number, Vertex());
    visitor_.insertVertices(position, number);
    return position;
}

/// Insert an additional edge.
///
/// \param vertexIndex0 Integer index of the first vertex in the edge.
/// \param vertexIndex1 Integer index of the second vertex in the edge.
/// \return Integer index of the newly inserted edge.
/// 
template<typename VISITOR>
inline size_t 
Graph<VISITOR>::insertEdge(
    const size_t vertexIndex0,
    const size_t vertexIndex1
) {
    assert(vertexIndex0 < numberOfVertices()); 
    assert(vertexIndex1 < numberOfVertices()); 
    
    if(multipleEdgesEnabled()) {
insertEdgeMark:
        edges_.push_back(Edge(vertexIndex0, vertexIndex1));
        size_t edgeIndex = edges_.size() - 1;
        insertAdjacenciesForEdge(edgeIndex);
        visitor_.insertEdge(edgeIndex);  
        return edgeIndex;
    }
    else {
        std::pair<bool, size_t> p = findEdge(vertexIndex0, vertexIndex1);
        if(p.first) { // edge already exists
            return p.second;
        }
        else {
            goto insertEdgeMark;
        }
    }
}

/// Erase a vertex and all edges connecting this vertex.
///
/// \param vertexIndex Integer index of the vertex to be erased.
/// 
template<typename VISITOR>
void 
Graph<VISITOR>::eraseVertex(
    const size_t vertexIndex
) {
    assert(vertexIndex < numberOfVertices()); 

    // erase all edges connected to the vertex
    while(vertices_[vertexIndex].size() != 0) {
        eraseEdge(vertices_[vertexIndex].begin()->edge());
    }

    if(vertexIndex == numberOfVertices()-1) { // if the last vertex is to be erased        
        vertices_.pop_back(); // erase vertex
        visitor_.eraseVertex(vertexIndex);
    }
    else { // if a vertex is to be erased which is not the last vertex
        // move last vertex to the free position:

        // collect indices of edges affected by the move
        size_t movingVertexIndex = numberOfVertices() - 1;
        std::set<size_t> affectedEdgeIndices;
        for(VertexIterator it = vertices_[movingVertexIndex].begin();
        it != vertices_[movingVertexIndex].end(); ++it) {
            affectedEdgeIndices.insert(it->edge());
        }
        
        // for all affected edges:
        for(std::set<size_t>::const_iterator it = affectedEdgeIndices.begin(); 
        it != affectedEdgeIndices.end(); ++it) { 
            // remove adjacencies
            eraseAdjacenciesForEdge(*it);

            // adapt vertex labels
            for(size_t j=0; j<2; ++j) {
                if(edges_[*it][j] == movingVertexIndex) {
                    edges_[*it][j] = vertexIndex;
                }
            }
            // if(!(edges_[*it].directedness()) && edges_[*it][0] > edges_[*it][1]) {
            if(edges_[*it][0] > edges_[*it][1]) {
                std::swap(edges_[*it][0], edges_[*it][1]);
            }
        }

        // move vertex
        vertices_[vertexIndex] = vertices_[movingVertexIndex]; // copy
        vertices_.pop_back(); // erase

        // insert adjacencies for edges of moved vertex
        for(std::set<size_t>::const_iterator it = affectedEdgeIndices.begin(); 
        it != affectedEdgeIndices.end(); ++it) { 
            insertAdjacenciesForEdge(*it);
        }

        visitor_.eraseVertex(vertexIndex);
        visitor_.relabelVertex(movingVertexIndex, vertexIndex);
    }
}

/// Erase an edge.
///
/// \param edgeIndex Integer index of the edge to be erased.
/// 
template<typename VISITOR>
inline void 
Graph<VISITOR>::eraseEdge(
    const size_t edgeIndex
) {
    assert(edgeIndex < numberOfEdges()); 

    eraseAdjacenciesForEdge(edgeIndex);
    if(edgeIndex == numberOfEdges() - 1) { // if the last edge is erased
        edges_.pop_back(); // delete
        visitor_.eraseEdge(edgeIndex);
    }
    else { 
        size_t movingEdgeIndex = numberOfEdges() - 1;
        eraseAdjacenciesForEdge(movingEdgeIndex);
        edges_[edgeIndex] = edges_[movingEdgeIndex]; // copy
        edges_.pop_back(); // delete
        insertAdjacenciesForEdge(edgeIndex);
        visitor_.eraseEdge(edgeIndex);
        visitor_.relabelEdge(movingEdgeIndex, edgeIndex);
    }
}

/// Get an iterator to the beginning of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
/// 
/// \sa verticesFromVertexEnd()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::VertexIterator 
Graph<VISITOR>::verticesFromVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].begin(); 
}

/// Get an iterator to the end of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
/// 
/// \sa verticesFromVertexBegin()
/// 
template<typename VISITOR>
inline typename Graph<VISITOR>::VertexIterator 
Graph<VISITOR>::verticesFromVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].end(); 
}

/// Get an iterator to the beginning of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
/// 
/// \sa verticesToVertexEnd()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::VertexIterator 
Graph<VISITOR>::verticesToVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].begin(); 
}

/// Get an iterator to the end of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
/// 
/// \sa verticesToVertexBegin()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::VertexIterator 
Graph<VISITOR>::verticesToVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].end(); 
}

/// Get an iterator to the beginning of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexEnd()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::EdgeIterator 
Graph<VISITOR>::edgesFromVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].begin(); 
}

/// Get an iterator to the end of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexBegin()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::EdgeIterator 
Graph<VISITOR>::edgesFromVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].end(); 
}

/// Get an iterator to the beginning of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexEnd()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::EdgeIterator 
Graph<VISITOR>::edgesToVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].begin(); 
}

/// Get an iterator to the end of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexBegin()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::EdgeIterator 
Graph<VISITOR>::edgesToVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].end(); 
}

/// Get an iterator to the beginning of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::AdjacencyIterator 
Graph<VISITOR>::adjacenciesFromVertexBegin(
    const size_t vertex
) const {
    return vertices_[vertex].begin();
}

/// Get an iterator to the end of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::AdjacencyIterator 
Graph<VISITOR>::adjacenciesFromVertexEnd(
    const size_t vertex
) const {
    return vertices_[vertex].end();
}

/// Get an iterator to the beginning of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::AdjacencyIterator 
Graph<VISITOR>::adjacenciesToVertexBegin(
    const size_t vertex
) const {
    return vertices_[vertex].begin();
}

/// Get an iterator to the end of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<typename VISITOR>
inline typename Graph<VISITOR>::AdjacencyIterator 
Graph<VISITOR>::adjacenciesToVertexEnd(
    const size_t vertex
) const {
    return vertices_[vertex].end();
}

/// Reserve memory for at least the given total number of vertices.
///
/// \param number Total number of vertices.
///
template<typename VISITOR>
inline void 
Graph<VISITOR>::reserveVertices(
    const size_t number
) {
    vertices_.reserve(number);
}

/// Reserve memory for at least the given total number of edges.
///
/// \param number Total number of edges.
///
template<typename VISITOR>
inline void 
Graph<VISITOR>::reserveEdges(
    const size_t number
) {
    edges_.reserve(number);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline const Adjacency<>&
Graph<VISITOR>::adjacencyFromVertex(
    const size_t vertex, 
    const size_t j
) const {
    return vertices_[vertex][j];
}

/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline const Adjacency<>&
Graph<VISITOR>::adjacencyToVertex(
    const size_t vertex, 
    const size_t j
) const {
    return vertices_[vertex][j];
}

/// Search for an edge (in logarithmic time).
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \return if an edge from vertex0 to vertex1 exists, pair.first is true 
///     and pair.second is the index of such an edge. if no edge from vertex0
///     to vertex1 exists, pair.first is false and pair.second is undefined.
///
template<typename VISITOR>
inline std::pair<bool, size_t>
Graph<VISITOR>::findEdge(
    const size_t vertex0, 
    const size_t vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());

    size_t v0 = vertex0;
    size_t v1 = vertex1;
    if(numberOfEdgesFromVertex(vertex1) < numberOfEdgesFromVertex(vertex0)) {
        v0 = vertex1;
        v1 = vertex0;
    }
    VertexIterator it = std::lower_bound(
        verticesFromVertexBegin(v0),
        verticesFromVertexEnd(v0),
        v1
    ); // binary search
    if(it != verticesFromVertexEnd(v0) && *it == v1) {
        // access the corresponding edge in constant time
        const size_t j = std::distance(verticesFromVertexBegin(v0), it);
        return std::make_pair(true, edgeFromVertex(v0, j));
    }
    else {
        return std::make_pair(false, 0);
    }
}

/// Indicate if multiple edges are enabled.
///
/// \return true if multiple edges are enabled, false otherwise.
///
template<typename VISITOR>
inline bool
Graph<VISITOR>::multipleEdgesEnabled() const {
    return multipleEdgesEnabled_;
}

/// Indicate if multiple edges are enabled.
///
/// Enable multiple edges like this: graph.multipleEdgesEnabled() = true;
///
/// \return reference the a Boolean flag.
///
template<typename VISITOR>
inline bool&
Graph<VISITOR>::multipleEdgesEnabled() {
    return multipleEdgesEnabled_;
}

template<typename VISITOR>
inline void 
Graph<VISITOR>::insertAdjacenciesForEdge(
    const size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const size_t vertexIndex0 = edge[0];
    const size_t vertexIndex1 = edge[1];
    vertices_[vertexIndex0].insert(
        Adjacency(vertexIndex1, edgeIndex)
    );
    if(vertexIndex1 != vertexIndex0) {
        vertices_[vertexIndex1].insert(
            Adjacency(vertexIndex0, edgeIndex)
        );
    }
}

template<typename VISITOR>
inline void 
Graph<VISITOR>::eraseAdjacenciesForEdge(
    const size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const size_t vertexIndex0 = edge[0];
    const size_t vertexIndex1 = edge[1];
    Vertex& vertex0 = vertices_[vertexIndex0];
    Vertex& vertex1 = vertices_[vertexIndex1];

    Adjacency adj(vertexIndex1, edgeIndex);
    RandomAccessSet<Adjacency>::iterator it = vertex0.find(adj);
    assert(it != vertex0.end()); 
    vertex0.erase(it);
    
    if(vertexIndex1 != vertexIndex0) { // if not a self-edge
        adj.vertex() = vertexIndex0;
        it = vertex1.find(adj);
        assert(it != vertex1.end()); 
        vertex1.erase(it);
    }
}

// implementation of Digraph

/// Construct a directed graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline 
Digraph<VISITOR>::Digraph(
    const Visitor& visitor
)
:   vertices_(),
    edges_(),
    multipleEdgesEnabled_(false),
    visitor_(visitor)
{}

/// Construct a directed graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline 
Digraph<VISITOR>::Digraph(
    const size_t numberOfVertices,
    const Visitor& visitor
)
:   vertices_(numberOfVertices),
    edges_(),
    multipleEdgesEnabled_(false),
    visitor_(visitor)
{
    visitor_.insertVertices(0, numberOfVertices);
}

/// Get the number of vertices.
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::numberOfVertices() const { 
    return vertices_.size(); 
}

/// Get the number of edges.
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::numberOfEdges() const { 
    return edges_.size(); 
}

/// Get the number of edges that originate from a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeFromVertex()
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::numberOfEdgesFromVertex(
    const size_t vertex
) const { 
    return vertices_[vertex].to_.size();
}

/// Get the number of edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::numberOfEdgesToVertex(
    const size_t vertex
) const { 
    return vertices_[vertex].from_.size();
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<typename VISITOR>
inline size_t
Digraph<VISITOR>::vertexOfEdge(
    const size_t edge,
    const size_t j
) const {
    assert(j < 2);

    return edges_[edge][j];
}

/// Get the integer index of an edge that originates from a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline size_t
Digraph<VISITOR>::edgeFromVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex].to_[j].edge();
}

/// Get the integer index of an edge that is incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesToVertex()
///
template<typename VISITOR>
inline size_t
Digraph<VISITOR>::edgeToVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex].from_[j].edge();
}

/// Get the integer index of a vertex reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline size_t
Digraph<VISITOR>::vertexFromVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex].to_[j].vertex();
}

/// Get the integer index of a vertex from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline size_t
Digraph<VISITOR>::vertexToVertex(
    const size_t vertex,
    const size_t j
) const {
    return vertices_[vertex].from_[j].vertex();
}

/// Insert an additional vertex.
///
/// \return Integer index of the newly inserted vertex.
///
/// \sa insertVertices()
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::insertVertex() {
    vertices_.push_back(Vertex());
    visitor_.insertVertex(vertices_.size() - 1);
    return vertices_.size() - 1;
}

/// Insert additional vertices.
///
/// \param number Number of new vertices to be inserted.
/// \return Integer index of the first newly inserted vertex.
///
/// \sa insertVertex()
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::insertVertices(
    const size_t number
) {
    size_t position = vertices_.size();
    vertices_.insert(vertices_.end(), number, Vertex());
    visitor_.insertVertices(position, number);
    return position;
}

/// Insert an additional edge.
///
/// \param vertexIndex0 Integer index of the first vertex (source) of the edge.
/// \param vertexIndex1 Integer index of the second vertex (target) of the edge.
/// \return Integer index of the newly inserted edge.
///
template<typename VISITOR>
inline size_t 
Digraph<VISITOR>::insertEdge(
    const size_t vertexIndex0,
    const size_t vertexIndex1
) {
    assert(vertexIndex0 < numberOfVertices()); 
    assert(vertexIndex1 < numberOfVertices()); 
    
    if(multipleEdgesEnabled()) {
insertEdgeMark:
        edges_.push_back(Edge(vertexIndex0, vertexIndex1));
        size_t edgeIndex = edges_.size() - 1;
        insertAdjacenciesForEdge(edgeIndex);
        visitor_.insertEdge(edgeIndex);  
        return edgeIndex;
    }
    else {
        std::pair<bool, size_t> p = findEdge(vertexIndex0, vertexIndex1);
        if(p.first) { // edge already exists
            return p.second;
        }
        else {
            goto insertEdgeMark;
        }
    }
}

/// Erase a vertex and all edges connecting this vertex.
///
/// \param vertexIndex Integer index of the vertex to be erased.
///
template<typename VISITOR>
void 
Digraph<VISITOR>::eraseVertex(
    const size_t vertexIndex
) {
    assert(vertexIndex < numberOfVertices()); 

    // erase all edges connected to the vertex
    while(vertices_[vertexIndex].from_.size() != 0) {
        eraseEdge(vertices_[vertexIndex].from_.begin()->edge());
    }
    while(vertices_[vertexIndex].to_.size() != 0) {
        eraseEdge(vertices_[vertexIndex].to_.begin()->edge());
    }

    if(vertexIndex == numberOfVertices() - 1) { // if the last vertex is to be erased        
        vertices_.pop_back(); // erase vertex
        visitor_.eraseVertex(vertexIndex);
    }
    else { // if a vertex is to be erased which is not the last vertex
        // move last vertex to the free position:

        // collect indices of edges affected by the move
        size_t movingVertexIndex = numberOfVertices() - 1;
        std::set<size_t> affectedEdgeIndices;
        for(VertexIterator it = vertices_[movingVertexIndex].from_.begin();
        it != vertices_[movingVertexIndex].from_.end(); ++it) {
            affectedEdgeIndices.insert(it->edge());
        }
        for(VertexIterator it = vertices_[movingVertexIndex].to_.begin();
        it != vertices_[movingVertexIndex].to_.end(); ++it) {
            affectedEdgeIndices.insert(it->edge());
        }
        
        // for all affected edges:
        for(std::set<size_t>::const_iterator it = affectedEdgeIndices.begin(); 
        it != affectedEdgeIndices.end(); ++it) { 
            // remove adjacencies
            eraseAdjacenciesForEdge(*it);

            // adapt vertex labels
            for(size_t j=0; j<2; ++j) {
                if(edges_[*it][j] == movingVertexIndex) {
                    edges_[*it][j] = vertexIndex;
                }
            }
        }

        // move vertex
        vertices_[vertexIndex] = vertices_[movingVertexIndex]; // copy
        vertices_.pop_back(); // erase

        // insert adjacencies for edges of moved vertex
        for(std::set<size_t>::const_iterator it = affectedEdgeIndices.begin(); 
        it != affectedEdgeIndices.end(); ++it) { 
            insertAdjacenciesForEdge(*it);
        }

        visitor_.eraseVertex(vertexIndex);
        visitor_.relabelVertex(movingVertexIndex, vertexIndex);
    }
}

/// Erase an edge.
///
/// \param edgeIndex Integer index of the edge to be erased.
///
template<typename VISITOR>
inline void 
Digraph<VISITOR>::eraseEdge(
    const size_t edgeIndex
) {
    assert(edgeIndex < numberOfEdges()); 

    eraseAdjacenciesForEdge(edgeIndex);
    if(edgeIndex == numberOfEdges() - 1) { // if the last edge is erased
        edges_.pop_back(); // delete
        visitor_.eraseEdge(edgeIndex);
    }
    else { 
        size_t movingEdgeIndex = numberOfEdges() - 1;
        eraseAdjacenciesForEdge(movingEdgeIndex);
        edges_[edgeIndex] = edges_[movingEdgeIndex]; // copy
        edges_.pop_back(); // delete
        insertAdjacenciesForEdge(edgeIndex);
        visitor_.eraseEdge(edgeIndex);
        visitor_.relabelEdge(movingEdgeIndex, edgeIndex);
    }
}

/// Get an iterator to the beginning of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexEnd()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::VertexIterator 
Digraph<VISITOR>::verticesFromVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].to_.begin(); 
}

/// Get an iterator to the end of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexBegin()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::VertexIterator 
Digraph<VISITOR>::verticesFromVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].to_.end(); 
}

/// Get an iterator to the beginning of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
/// 
/// \sa verticesToVertexEnd()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::VertexIterator 
Digraph<VISITOR>::verticesToVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].from_.begin(); 
}

/// Get an iterator to the end of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexBegin()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::VertexIterator 
Digraph<VISITOR>::verticesToVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].from_.end(); 
}

/// Get an iterator to the beginning of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexEnd()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::EdgeIterator 
Digraph<VISITOR>::edgesFromVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].to_.begin(); 
}

/// Get an iterator to the end of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexBegin()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::EdgeIterator 
Digraph<VISITOR>::edgesFromVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].to_.end(); 
}

/// Get an iterator to the beginning of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexEnd()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::EdgeIterator 
Digraph<VISITOR>::edgesToVertexBegin(
    const size_t vertex
) const { 
    return vertices_[vertex].from_.begin(); 
}

/// Get an iterator to the end of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexBegin()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::EdgeIterator 
Digraph<VISITOR>::edgesToVertexEnd(
    const size_t vertex
) const { 
    return vertices_[vertex].from_.end(); 
}

/// Get an iterator to the beginning of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::AdjacencyIterator 
Digraph<VISITOR>::adjacenciesFromVertexBegin(
    const size_t vertex
) const {
    return vertices_[vertex].to_.begin();
}

/// Get an iterator to the end of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::AdjacencyIterator 
Digraph<VISITOR>::adjacenciesFromVertexEnd(
    const size_t vertex
) const {
    return vertices_[vertex].to_.end();
}

/// Get an iterator to the beginning of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::AdjacencyIterator 
Digraph<VISITOR>::adjacenciesToVertexBegin(
    const size_t vertex
) const {
    return vertices_[vertex].from_.begin();
}

/// Get an iterator to the end of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<typename VISITOR>
inline typename Digraph<VISITOR>::AdjacencyIterator 
Digraph<VISITOR>::adjacenciesToVertexEnd(
    const size_t vertex
) const {
    return vertices_[vertex].from_.end();
}

/// Reserve memory for at least the given total number of vertices.
///
/// \param number Total number of vertices.
///
template<typename VISITOR>
inline void 
Digraph<VISITOR>::reserveVertices(
    const size_t number
) {
    vertices_.reserve(number);
}

/// Reserve memory for at least the given total number of edges.
///
/// \param number Total number of edges.
///
template<typename VISITOR>
inline void 
Digraph<VISITOR>::reserveEdges(
    const size_t number
) {
    edges_.reserve(number);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline const Adjacency<>&
Digraph<VISITOR>::adjacencyFromVertex(
    const size_t vertex, 
    const size_t j
) const {
    return vertices_[vertex].to_[j];
}

/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline const Adjacency<>&
Digraph<VISITOR>::adjacencyToVertex(
    const size_t vertex, 
    const size_t j
) const {
    return vertices_[vertex].from_[j];
}

/// Search for an edge (in logarithmic time).
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \return if an edge from vertex0 to vertex1 exists, pair.first is true 
///     and pair.second is the index of such an edge. if no edge from vertex0
///     to vertex1 exists, pair.first is false and pair.second is undefined.
///
template<typename VISITOR>
inline std::pair<bool, size_t>
Digraph<VISITOR>::findEdge(
    const size_t vertex0, 
    const size_t vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());

    VertexIterator it = std::lower_bound(
        verticesFromVertexBegin(vertex0),
        verticesFromVertexEnd(vertex0),
        vertex1
    ); // binary search
    if(it != verticesFromVertexEnd(vertex0) && *it == vertex1) {
        // access the corresponding edge in constant time
        const size_t j = std::distance(verticesFromVertexBegin(vertex0), it);
        return std::make_pair(true, edgeFromVertex(vertex0, j));
    }
    else {
        return std::make_pair(false, 0);
    }
}

/// Indicate if multiple edges are enabled.
///
/// \return true if multiple edges are enabled, false otherwise.
///
template<typename VISITOR>
inline bool
Digraph<VISITOR>::multipleEdgesEnabled() const {
    return multipleEdgesEnabled_;
}

/// Indicate if multiple edges are enabled.
///
/// Enable multiple edges like this: graph.multipleEdgesEnabled() = true;
///
/// \return reference the a Boolean flag.
///
template<typename VISITOR>
inline bool&
Digraph<VISITOR>::multipleEdgesEnabled() {
    return multipleEdgesEnabled_;
}

template<typename VISITOR>
inline void 
Digraph<VISITOR>::insertAdjacenciesForEdge(
    const size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const size_t vertexIndex0 = edge[0];
    const size_t vertexIndex1 = edge[1];
    vertices_[vertexIndex0].to_.insert(
        Adjacency(vertexIndex1, edgeIndex)
    );
    vertices_[vertexIndex1].from_.insert(
        Adjacency(vertexIndex0, edgeIndex)
    );
}

template<typename VISITOR>
inline void 
Digraph<VISITOR>::eraseAdjacenciesForEdge(
    const size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const size_t vertexIndex0 = edge[0];
    const size_t vertexIndex1 = edge[1];
    Vertex& vertex0 = vertices_[vertexIndex0];
    Vertex& vertex1 = vertices_[vertexIndex1];

    Adjacency adj(vertexIndex1, edgeIndex);
    RandomAccessSet<Adjacency>::iterator it = vertex0.to_.find(adj);
    assert(it != vertex0.to_.end()); 
    vertex0.to_.erase(it);
    
    adj.vertex() = vertexIndex0;
    it = vertex1.from_.find(adj);
    assert(it != vertex1.from_.end()); 
    vertex1.from_.erase(it);
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_HXX
