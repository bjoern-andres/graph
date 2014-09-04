/// \mainpage
/// Graphs and Graph Algorithms in C++.
///
/// Copyright (c) 2013 by Bjoern Andres, http://www.andres.sc
///
/// This software was developed by Bjoern Andres, Duligur Ibeling,
/// and Mark Matten.
/// Enquiries shall be directed to bjoern@andres.sc.
///
/// \section section_abstract Short Description
/// This set of header files implements directed and undirected graphs as 
/// adjacency lists with constant-time access. Vertices and edges are always 
/// indexed by contiguous integers. This indexing simplifies the implementation 
/// of algorithms for static graphs. In dynamic settings where vertices and 
/// edges are removed from a graph, indices of vertices and edges can change. 
/// These changes can be followed, if necessary, by means of a visitor.
///
/// \section section_features Features
/// - Directed and undirected graphs, implemented as adjacency lists
///   with constant-time access
/// - Access to vertices and edges by contiguous integer indices
/// - Access to neighboring vertices and incident edges by STL-compliant 
///   random access iterators
/// - Insertion and removal of vertices and edges 
/// - Multiple edges, which are disabled by default, can be enabled
/// - Visitors that follow changes of vertex and edge indices 
/// - Graph algorithms
///   - Connected components 
///     - by breadth-first search
///     - by disjoint sets 
///     .
///   - Shortest paths in weighted and unweighted graphs, as sequences of edges 
///     or vertices
///     - Single source shortest path (SSSP)
///     - Single pair shortest path (SPSP)
///     .
///   - Maximum s-t-flow
///     - Push-relabel algorithm with FIFO vertex selection rule.
///     - Edmonds-Karp algorithm
///     .
///   - Minimum multicuts by interger programming, using Cplex or Gurobi 
///   .
/// 
/// \section section_license License
/// Copyright (c) 2013 by Bjoern Andres.
/// 
/// This software was developed by Bjoern Andres, Duligur Ibeling,
/// and Mark Matten.
/// Enquiries shall be directed to bjoern@andres.sc.
///
/// Redistribution and use in source and binary forms, with or without 
/// modification, are permitted provided that the following conditions are met:
/// - Redistributions of source code must retain the above copyright notice,
///   this list of conditions and the following disclaimer.
/// - Redistributions in binary form must reproduce the above copyright notice, 
///   this list of conditions and the following disclaimer in the documentation
///   and/or other materials provided with the distribution.
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
#include <cstddef>
#include <cmath>
#include <iterator> // std::random_access_iterator
#include <vector>
#include <set> 
#include <iostream>
#include <utility> // std::pair
#include <algorithm> // std::fill

#include "andres/random-access-set.hxx"

/// The public API.
namespace andres {

/// Graphs and Algorithms that Operate on Graphs.
namespace graph {

/// The adjacency of a vertex consists of a vertex and a connecting edge.
template<class T = std::size_t>
class Adjacency {
public:
    typedef T Value;

    Adjacency(const Value = T(), const Value = T());
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
template<class T = std::size_t>
struct DefaultSubgraphMask {
    typedef T Value;

    bool vertex(const Value v) const { return true; }
    bool edge(const Value e) const { return true; }
};

/// Return 1 for every edge.
template<class T>
struct UnitEdgeWeightIterator {
    typedef std::ptrdiff_t difference_type;
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::random_access_iterator_tag iterator_category;

    value_type operator[](const std::size_t) const
        { return 1; }
};

// \cond SUPPRESS_DOXYGEN
namespace graph_detail {

template<bool DIRECTED, class T = std::size_t>
class Edge {
public:
    typedef T Value;

    Edge(const Value, const Value);
    Value operator[](const std::size_t) const;
    Value& operator[](const std::size_t);

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
    #ifdef _MSC_VER
    difference_type operator-(const IteratorHelper<T>&) const;
    #endif

    // access
    std::size_t operator*() const;
    std::size_t operator[](const std::size_t j) const;
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
    void insertVertex(const std::size_t a) const {}
    void insertVertices(const std::size_t a, const std::size_t n) const {}
    void eraseVertex(const std::size_t a) const {}
    void relabelVertex(const std::size_t a, const std::size_t b) const {}
    void insertEdge(const std::size_t a) const {}
    void eraseEdge(const std::size_t a) const {}
    void relabelEdge(const std::size_t a, const std::size_t b) const {}
};

/// Visitors can be used to follow the indices of vertices and edges.
///
/// These indices change due to the insetion and removal of vertices and edges.
///
struct VerboseGraphVisitor {
    VerboseGraphVisitor() {}
    void insertVertex(const std::size_t a) const
        { std::cout << "inserting vertex " << a << std::endl; }
    void insertVertices(const std::size_t a, const std::size_t n) const
        { std::cout << "inserting " << n << " vertices, starting from index " << a << std::endl; }
    void eraseVertex(const std::size_t a) const
        { std::cout << "removing vertex " << a << std::endl; }
    void relabelVertex(const std::size_t a, const std::size_t b) const
        { std::cout << "relabeling vertex " << a << ". new label is " << b << std::endl; }
    void insertEdge(const std::size_t a) const
        { std::cout << "inserting edge " << a << std::endl; }
    void eraseEdge(const std::size_t a) const
        { std::cout << "removing edge " << a << std::endl; }
    void relabelEdge(const std::size_t a, const std::size_t b) const
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
    Graph(const std::size_t, const Visitor& = Visitor());
    void assign(const Visitor& = Visitor());
    void assign(const std::size_t, const Visitor& = Visitor());
    void reserveVertices(const std::size_t);
    void reserveEdges(const std::size_t);

    // iterator access (compatible with Digraph)
    VertexIterator verticesFromVertexBegin(const std::size_t) const;
    VertexIterator verticesFromVertexEnd(const std::size_t) const;
    VertexIterator verticesToVertexBegin(const std::size_t) const;
    VertexIterator verticesToVertexEnd(const std::size_t) const;
    EdgeIterator edgesFromVertexBegin(const std::size_t) const;
    EdgeIterator edgesFromVertexEnd(const std::size_t) const;
    EdgeIterator edgesToVertexBegin(const std::size_t) const;
    EdgeIterator edgesToVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexEnd(const std::size_t) const;

    // access (compatible with Digraph)
    std::size_t numberOfVertices() const;
    std::size_t numberOfEdges() const;
    std::size_t numberOfEdgesFromVertex(const std::size_t) const;
    std::size_t numberOfEdgesToVertex(const std::size_t) const;
    std::size_t vertexOfEdge(const std::size_t, const std::size_t) const;
    std::size_t edgeFromVertex(const std::size_t, const std::size_t) const;
    std::size_t edgeToVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexFromVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexToVertex(const std::size_t, const std::size_t) const;
    const Adjacency<>& adjacencyFromVertex(const std::size_t, const std::size_t) const;
    const Adjacency<>& adjacencyToVertex(const std::size_t, const std::size_t) const;
    std::pair<bool, std::size_t> findEdge(const std::size_t, const std::size_t) const;
    bool multipleEdgesEnabled() const;

    // manipulation
    std::size_t insertVertex();
    std::size_t insertVertices(const std::size_t);
    std::size_t insertEdge(const std::size_t, const std::size_t);
    void eraseVertex(const std::size_t);
    void eraseEdge(const std::size_t);
    bool& multipleEdgesEnabled();

private:
    typedef Adjacency<> AdjacencyType;
    typedef graph_detail::Adjacencies Vertex;
    typedef graph_detail::Edge<false> Edge;

    void insertAdjacenciesForEdge(const std::size_t);
    void eraseAdjacenciesForEdge(const std::size_t);

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
    Digraph(const std::size_t, const Visitor& = Visitor());
    void assign(const Visitor& = Visitor());
    void assign(const std::size_t, const Visitor& = Visitor());
    void reserveVertices(const std::size_t);
    void reserveEdges(const std::size_t);

    // iterator access 
    VertexIterator verticesFromVertexBegin(const std::size_t) const;
    VertexIterator verticesFromVertexEnd(const std::size_t) const;
    VertexIterator verticesToVertexBegin(const std::size_t) const;
    VertexIterator verticesToVertexEnd(const std::size_t) const;
    EdgeIterator edgesFromVertexBegin(const std::size_t) const;
    EdgeIterator edgesFromVertexEnd(const std::size_t) const;
    EdgeIterator edgesToVertexBegin(const std::size_t) const;
    EdgeIterator edgesToVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexEnd(const std::size_t) const;

    // access 
    std::size_t numberOfVertices() const;
    std::size_t numberOfEdges() const;
    std::size_t numberOfEdgesFromVertex(const std::size_t) const;
    std::size_t numberOfEdgesToVertex(const std::size_t) const;
    std::size_t vertexOfEdge(const std::size_t, const std::size_t) const;
    std::size_t edgeFromVertex(const std::size_t, const std::size_t) const;
    std::size_t edgeToVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexFromVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexToVertex(const std::size_t, const std::size_t) const;
    const Adjacency<>& adjacencyFromVertex(const std::size_t, const std::size_t) const;
    const Adjacency<>& adjacencyToVertex(const std::size_t, const std::size_t) const;
    std::pair<bool, std::size_t> findEdge(const std::size_t, const std::size_t) const;
    bool multipleEdgesEnabled() const;

    // manipulation
    std::size_t insertVertex();
    std::size_t insertVertices(const std::size_t);
    std::size_t insertEdge(const std::size_t, const std::size_t);
    void eraseVertex(const std::size_t);
    void eraseEdge(const std::size_t);
    bool& multipleEdgesEnabled();

private:
    typedef Adjacency<> AdjacencyType;
    typedef graph_detail::Adjacencies Adjacencies;
    struct Vertex {
        Vertex() 
            : from_(), to_() 
            {}
        Adjacencies from_;
        Adjacencies to_;
    };
    typedef graph_detail::Edge<true> Edge;

    void insertAdjacenciesForEdge(const std::size_t);
    void eraseAdjacenciesForEdge(const std::size_t);

    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    bool multipleEdgesEnabled_;
    Visitor visitor_;
};

/// Complete graph.
template<typename VISITOR = IdleGraphVisitor>
class CompleteGraph {
public:
    typedef VISITOR Visitor;

    // \cond SUPPRESS_DOXYGEN
    class AdjacencyIterator {
    public:
        typedef CompleteGraph<Visitor> GraphType;
        typedef typename std::random_access_iterator_tag iterator_category;
        typedef Adjacency<> value_type;
        typedef typename std::ptrdiff_t difference_type;

        AdjacencyIterator();
        AdjacencyIterator(const GraphType&);
        AdjacencyIterator(const GraphType&, const std::size_t);
        AdjacencyIterator(const GraphType&, const std::size_t, const std::size_t);

        // increment and decrement
        AdjacencyIterator& operator+=(const difference_type);
        AdjacencyIterator& operator-=(const difference_type);
        AdjacencyIterator& operator++(); // prefix
        AdjacencyIterator& operator--(); // prefix
        AdjacencyIterator operator++(int); // postfix
        AdjacencyIterator operator--(int); // postfix
        AdjacencyIterator operator+(const difference_type) const;
        AdjacencyIterator operator-(const difference_type) const;

        // comparison
        bool operator==(const AdjacencyIterator&) const;
        bool operator!=(const AdjacencyIterator&) const;
        bool operator<(const AdjacencyIterator&) const;
        bool operator<=(const AdjacencyIterator&) const;
        bool operator>(const AdjacencyIterator&) const;
        bool operator>=(const AdjacencyIterator&) const;

        // access
        const Adjacency<>& operator*();
        const Adjacency<>* operator->();
        const Adjacency<>& operator[](const std::size_t);

    protected:
        const GraphType* graph_;
        std::size_t vertex_;
        std::size_t adjacencyIndex_;
        Adjacency<> adjacency_;
    };

    class VertexIterator
    :   public AdjacencyIterator {
    public:
        typedef CompleteGraph<Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef typename Base::iterator_category iterator_category;
        typedef typename Base::value_type value_type;
        typedef typename Base::difference_type difference_type;

        VertexIterator();
        VertexIterator(const VertexIterator&);
        VertexIterator(const AdjacencyIterator&);
        VertexIterator(const GraphType&);
        VertexIterator(const GraphType&, const std::size_t);
        VertexIterator(const GraphType&, const std::size_t, const std::size_t);

        // access
        std::size_t operator*() const;
        std::size_t operator[](const std::size_t) const;
    };

    class EdgeIterator
    :   public AdjacencyIterator {
    public:
        typedef CompleteGraph<Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef typename Base::iterator_category iterator_category;
        typedef typename Base::value_type value_type;
        typedef typename Base::difference_type difference_type;

        EdgeIterator();
        EdgeIterator(const EdgeIterator&);
        EdgeIterator(const AdjacencyIterator&);
        EdgeIterator(const GraphType&);
        EdgeIterator(const GraphType&, const std::size_t);
        EdgeIterator(const GraphType&, const std::size_t, const std::size_t);

        // access
        std::size_t operator*() const;
        std::size_t operator[](const std::size_t) const;
    };
    // \endcond

    // construction
    CompleteGraph(const Visitor& = Visitor());
    CompleteGraph(const std::size_t, const Visitor& = Visitor());
    void assign(const Visitor& = Visitor());
    void assign(const std::size_t, const Visitor& = Visitor());

    // iterator access (compatible with Digraph)
    VertexIterator verticesFromVertexBegin(const std::size_t) const;
    VertexIterator verticesFromVertexEnd(const std::size_t) const;
    VertexIterator verticesToVertexBegin(const std::size_t) const;
    VertexIterator verticesToVertexEnd(const std::size_t) const;
    EdgeIterator edgesFromVertexBegin(const std::size_t) const;
    EdgeIterator edgesFromVertexEnd(const std::size_t) const;
    EdgeIterator edgesToVertexBegin(const std::size_t) const;
    EdgeIterator edgesToVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexEnd(const std::size_t) const;

    // access (compatible with Digraph)
    std::size_t numberOfVertices() const;
    std::size_t numberOfEdges() const;
    std::size_t numberOfEdgesFromVertex(const std::size_t) const;
    std::size_t numberOfEdgesToVertex(const std::size_t) const;
    std::size_t vertexOfEdge(const std::size_t, const std::size_t) const;
    std::size_t edgeFromVertex(const std::size_t, const std::size_t) const;
    std::size_t edgeToVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexFromVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexToVertex(const std::size_t, const std::size_t) const;
    Adjacency<> adjacencyFromVertex(const std::size_t, const std::size_t) const;
    Adjacency<> adjacencyToVertex(const std::size_t, const std::size_t) const;
    std::pair<bool, std::size_t> findEdge(const std::size_t, const std::size_t) const;
    bool multipleEdgesEnabled() const;

private:
    std::size_t edgeOfStrictlyIncreasingPairOfVertices(const std::size_t, const std::size_t) const;

    std::size_t numberOfVertices_;
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
    const std::size_t j
) const {
    assert(j < 2);

    return vertexIndices_[j];
}

template<bool DIRECTED, class T>
inline typename Edge<DIRECTED, T>::Value&
Edge<DIRECTED, T>::operator[](
    const std::size_t j
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
inline std::size_t
IteratorHelper<T>::operator*() const { 
    if(T) { // evaluated at compile time
        return Base::operator*().vertex(); 
    }
    else {
        return Base::operator*().edge(); 
    }
}

template<bool T>
inline std::size_t
IteratorHelper<T>::operator[](
    const std::size_t j
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

#ifdef _MSC_VER
template<bool T>
inline typename IteratorHelper<T>::difference_type
IteratorHelper<T>::operator-(
    const IteratorHelper<T>& other
) const {
    return Base::operator-(other);
}
#endif

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
    const std::size_t numberOfVertices,
    const Visitor& visitor
)
:   vertices_(numberOfVertices),
    edges_(),
    multipleEdgesEnabled_(false),
    visitor_(visitor)
{
    visitor_.insertVertices(0, numberOfVertices);
}

/// Clear an undirected graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
Graph<VISITOR>::assign(
    const Visitor& visitor
) {
    vertices_.clear();
    edges_.clear();
    multipleEdgesEnabled_ = false;
    visitor_ = visitor;
}

/// Clear an undirected graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
Graph<VISITOR>::assign(
    const std::size_t numberOfVertices,
    const Visitor& visitor
) {
    vertices_.resize(numberOfVertices);
    std::fill(vertices_.begin(), vertices_.end(), Vertex());
    edges_.clear();
    multipleEdgesEnabled_ = false;
    visitor_ = visitor;
    visitor_.insertVertices(0, numberOfVertices);
}
    
/// Get the number of vertices.
///
template<typename VISITOR>
inline std::size_t
Graph<VISITOR>::numberOfVertices() const { 
    return vertices_.size(); 
}

/// Get the number of edges.
///
template<typename VISITOR>
inline std::size_t
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
inline std::size_t
Graph<VISITOR>::numberOfEdgesFromVertex(
    const std::size_t vertex
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
inline std::size_t
Graph<VISITOR>::numberOfEdgesToVertex(
    const std::size_t vertex
) const { 
    return vertices_[vertex].size();
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<typename VISITOR>
inline std::size_t
Graph<VISITOR>::vertexOfEdge(
    const std::size_t edge,
    const std::size_t j
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
inline std::size_t
Graph<VISITOR>::edgeFromVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
Graph<VISITOR>::edgeToVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
Graph<VISITOR>::vertexFromVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
Graph<VISITOR>::vertexToVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
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
inline std::size_t
Graph<VISITOR>::insertVertices(
    const std::size_t number
) {
    std::size_t position = vertices_.size();
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
inline std::size_t
Graph<VISITOR>::insertEdge(
    const std::size_t vertexIndex0,
    const std::size_t vertexIndex1
) {
    assert(vertexIndex0 < numberOfVertices()); 
    assert(vertexIndex1 < numberOfVertices()); 
    
    if(multipleEdgesEnabled()) {
insertEdgeMark:
        edges_.push_back(Edge(vertexIndex0, vertexIndex1));
        std::size_t edgeIndex = edges_.size() - 1;
        insertAdjacenciesForEdge(edgeIndex);
        visitor_.insertEdge(edgeIndex);  
        return edgeIndex;
    }
    else {
        std::pair<bool, std::size_t> p = findEdge(vertexIndex0, vertexIndex1);
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
    const std::size_t vertexIndex
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
        std::size_t movingVertexIndex = numberOfVertices() - 1;
        std::set<std::size_t> affectedEdgeIndices;
        for(Vertex::const_iterator it = vertices_[movingVertexIndex].begin();
        it != vertices_[movingVertexIndex].end(); ++it) {
            affectedEdgeIndices.insert(it->edge());
        }
        
        // for all affected edges:
        for(std::set<std::size_t>::const_iterator it = affectedEdgeIndices.begin();
        it != affectedEdgeIndices.end(); ++it) { 
            // remove adjacencies
            eraseAdjacenciesForEdge(*it);

            // adapt vertex labels
            for(std::size_t j=0; j<2; ++j) {
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
        for(std::set<std::size_t>::const_iterator it = affectedEdgeIndices.begin();
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
    const std::size_t edgeIndex
) {
    assert(edgeIndex < numberOfEdges()); 

    eraseAdjacenciesForEdge(edgeIndex);
    if(edgeIndex == numberOfEdges() - 1) { // if the last edge is erased
        edges_.pop_back(); // delete
        visitor_.eraseEdge(edgeIndex);
    }
    else { 
        std::size_t movingEdgeIndex = numberOfEdges() - 1;
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t number
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
    const std::size_t number
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
    const std::size_t vertex,
    const std::size_t j
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
    const std::size_t vertex,
    const std::size_t j
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
inline std::pair<bool, std::size_t>
Graph<VISITOR>::findEdge(
    const std::size_t vertex0,
    const std::size_t vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());

    std::size_t v0 = vertex0;
    std::size_t v1 = vertex1;
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
        const std::size_t j = std::distance(verticesFromVertexBegin(v0), it);
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
    const std::size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const std::size_t vertexIndex0 = edge[0];
    const std::size_t vertexIndex1 = edge[1];
    vertices_[vertexIndex0].insert(
        AdjacencyType(vertexIndex1, edgeIndex)
    );
    if(vertexIndex1 != vertexIndex0) {
        vertices_[vertexIndex1].insert(
            AdjacencyType(vertexIndex0, edgeIndex)
        );
    }
}

template<typename VISITOR>
inline void 
Graph<VISITOR>::eraseAdjacenciesForEdge(
    const std::size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const std::size_t vertexIndex0 = edge[0];
    const std::size_t vertexIndex1 = edge[1];
    Vertex& vertex0 = vertices_[vertexIndex0];
    Vertex& vertex1 = vertices_[vertexIndex1];

    AdjacencyType adj(vertexIndex1, edgeIndex);
    RandomAccessSet<AdjacencyType>::iterator it = vertex0.find(adj);
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
    const std::size_t numberOfVertices,
    const Visitor& visitor
)
:   vertices_(numberOfVertices),
    edges_(),
    multipleEdgesEnabled_(false),
    visitor_(visitor)
{
    visitor_.insertVertices(0, numberOfVertices);
}

/// Clear a directed graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
Digraph<VISITOR>::assign(
    const Visitor& visitor
) {
    vertices_.clear();
    edges_.clear();
    multipleEdgesEnabled_ = false;
    visitor_ = visitor;
}

/// Clear a directed graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
Digraph<VISITOR>::assign(
    const std::size_t numberOfVertices,
    const Visitor& visitor
) {
    vertices_.resize(numberOfVertices);
    std::fill(vertices_.begin(), vertices_.end(), Vertex());
    edges_.clear();
    multipleEdgesEnabled_ = false;
    visitor_ = visitor;
    visitor_.insertVertices(0, numberOfVertices);
}

/// Get the number of vertices.
///
template<typename VISITOR>
inline std::size_t
Digraph<VISITOR>::numberOfVertices() const { 
    return vertices_.size(); 
}

/// Get the number of edges.
///
template<typename VISITOR>
inline std::size_t
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
inline std::size_t
Digraph<VISITOR>::numberOfEdgesFromVertex(
    const std::size_t vertex
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
inline std::size_t
Digraph<VISITOR>::numberOfEdgesToVertex(
    const std::size_t vertex
) const { 
    return vertices_[vertex].from_.size();
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<typename VISITOR>
inline std::size_t
Digraph<VISITOR>::vertexOfEdge(
    const std::size_t edge,
    const std::size_t j
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
inline std::size_t
Digraph<VISITOR>::edgeFromVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
Digraph<VISITOR>::edgeToVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
Digraph<VISITOR>::vertexFromVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
Digraph<VISITOR>::vertexToVertex(
    const std::size_t vertex,
    const std::size_t j
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
inline std::size_t
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
inline std::size_t
Digraph<VISITOR>::insertVertices(
    const std::size_t number
) {
    std::size_t position = vertices_.size();
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
inline std::size_t
Digraph<VISITOR>::insertEdge(
    const std::size_t vertexIndex0,
    const std::size_t vertexIndex1
) {
    assert(vertexIndex0 < numberOfVertices()); 
    assert(vertexIndex1 < numberOfVertices()); 
    
    if(multipleEdgesEnabled()) {
insertEdgeMark:
        edges_.push_back(Edge(vertexIndex0, vertexIndex1));
        std::size_t edgeIndex = edges_.size() - 1;
        insertAdjacenciesForEdge(edgeIndex);
        visitor_.insertEdge(edgeIndex);  
        return edgeIndex;
    }
    else {
        std::pair<bool, std::size_t> p = findEdge(vertexIndex0, vertexIndex1);
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
    const std::size_t vertexIndex
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
        std::size_t movingVertexIndex = numberOfVertices() - 1;
        std::set<std::size_t> affectedEdgeIndices;
        for(Adjacencies::const_iterator it = vertices_[movingVertexIndex].from_.begin();
        it != vertices_[movingVertexIndex].from_.end(); ++it) {
            affectedEdgeIndices.insert(it->edge());
        }
        for(Adjacencies::const_iterator it = vertices_[movingVertexIndex].to_.begin();
        it != vertices_[movingVertexIndex].to_.end(); ++it) {
            affectedEdgeIndices.insert(it->edge());
        }
        
        // for all affected edges:
        for(std::set<std::size_t>::const_iterator it = affectedEdgeIndices.begin();
        it != affectedEdgeIndices.end(); ++it) { 
            // remove adjacencies
            eraseAdjacenciesForEdge(*it);

            // adapt vertex labels
            for(std::size_t j=0; j<2; ++j) {
                if(edges_[*it][j] == movingVertexIndex) {
                    edges_[*it][j] = vertexIndex;
                }
            }
        }

        // move vertex
        vertices_[vertexIndex] = vertices_[movingVertexIndex]; // copy
        vertices_.pop_back(); // erase

        // insert adjacencies for edges of moved vertex
        for(std::set<std::size_t>::const_iterator it = affectedEdgeIndices.begin();
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
    const std::size_t edgeIndex
) {
    assert(edgeIndex < numberOfEdges()); 

    eraseAdjacenciesForEdge(edgeIndex);
    if(edgeIndex == numberOfEdges() - 1) { // if the last edge is erased
        edges_.pop_back(); // delete
        visitor_.eraseEdge(edgeIndex);
    }
    else { 
        std::size_t movingEdgeIndex = numberOfEdges() - 1;
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t vertex
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
    const std::size_t number
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
    const std::size_t number
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
    const std::size_t vertex,
    const std::size_t j
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
    const std::size_t vertex,
    const std::size_t j
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
inline std::pair<bool, std::size_t>
Digraph<VISITOR>::findEdge(
    const std::size_t vertex0,
    const std::size_t vertex1
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
        const std::size_t j = std::distance(verticesFromVertexBegin(vertex0), it);
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
    const std::size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const std::size_t vertexIndex0 = edge[0];
    const std::size_t vertexIndex1 = edge[1];
    vertices_[vertexIndex0].to_.insert(
        AdjacencyType(vertexIndex1, edgeIndex)
    );
    vertices_[vertexIndex1].from_.insert(
        AdjacencyType(vertexIndex0, edgeIndex)
    );
}

template<typename VISITOR>
inline void 
Digraph<VISITOR>::eraseAdjacenciesForEdge(
    const std::size_t edgeIndex
) {
    const Edge& edge = edges_[edgeIndex];
    const std::size_t vertexIndex0 = edge[0];
    const std::size_t vertexIndex1 = edge[1];
    Vertex& vertex0 = vertices_[vertexIndex0];
    Vertex& vertex1 = vertices_[vertexIndex1];

    AdjacencyType adj(vertexIndex1, edgeIndex);
    RandomAccessSet<AdjacencyType>::iterator it = vertex0.to_.find(adj);
    assert(it != vertex0.to_.end()); 
    vertex0.to_.erase(it);
    
    adj.vertex() = vertexIndex0;
    it = vertex1.from_.find(adj);
    assert(it != vertex1.from_.end()); 
    vertex1.from_.erase(it);
}

/// Construct a complete graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline
CompleteGraph<VISITOR>::CompleteGraph(
    const Visitor& visitor
)
:   numberOfVertices_(0),
    visitor_(visitor)
{}

/// Construct a complete graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline
CompleteGraph<VISITOR>::CompleteGraph(
    const std::size_t numberOfVertices,
    const Visitor& visitor
)
:   numberOfVertices_(numberOfVertices),
    visitor_(visitor)
{}

/// Clear a complete graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
CompleteGraph<VISITOR>::assign(
    const Visitor& visitor
) {
    numberOfVertices_ = 0;
    visitor_ = visitor;
}

/// Clear a complete graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
CompleteGraph<VISITOR>::assign(
    const std::size_t numberOfVertices,
    const Visitor& visitor
) {
    numberOfVertices_ = numberOfVertices;
    visitor_ = visitor;
}

/// Get an iterator to the beginning of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesFromVertexBegin(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesFromVertexEnd(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesToVertexBegin(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesToVertexEnd(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesFromVertexBegin(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesFromVertexEnd(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesToVertexBegin(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesToVertexEnd(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesFromVertexBegin(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesFromVertexEnd(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesToVertexBegin(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesToVertexEnd(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get the number of vertices.
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfVertices() const {
    return numberOfVertices_;
}

/// Get the number of edges.
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfEdges() const {
    return numberOfVertices() * (numberOfVertices() - 1) / 2;
}

/// Get the number of edges that originate from a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfEdgesFromVertex(
    const std::size_t vertex
) const {
    assert(vertex < numberOfVertices());
    return numberOfVertices() - 1;
}

/// Get the number of edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfEdgesToVertex(
    const std::size_t vertex
) const {
    assert(vertex < numberOfVertices());
    return numberOfVertices() - 1;
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::vertexOfEdge(
    const std::size_t edge,
    const std::size_t j
) const {
    assert(edge < numberOfEdges());
    assert(j < 2);
    const float p = static_cast<float>(numberOfVertices() * 2 - 1) / 2;
    const float q = static_cast<float>(edge) * 2;
    const std::size_t vertex0 = static_cast<std::size_t>(p - std::sqrt(p * p - q));
    if(j == 0) {
        return vertex0;
    }
    else {
        return edge + vertex0 * (vertex0 + 1) / 2 - (numberOfVertices() - 1) * vertex0 + 1;
    }
}

/// Get the integer index of an edge that originates from a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::edgeFromVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesFromVertex(vertex));
    if(j < vertex) {
        const std::size_t vertexAdjacent = j;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertexAdjacent, vertex);
        return edgeAdjacent;
    }
    else {
        const std::size_t vertexAdjacent = j + 1;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertex, vertexAdjacent);
        return edgeAdjacent;
    }
}

/// Get the integer index of an edge that is incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesToVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::edgeToVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return edgeFromVertex(vertex, j);
}

/// Get the integer index of a vertex reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::vertexFromVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesFromVertex(vertex));
    if(j < vertex) {
        return j;
    }
    else {
        return j + 1;
    }
}

/// Get the integer index of a vertex from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::vertexToVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return vertexFromVertex(vertex, j);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline Adjacency<>
CompleteGraph<VISITOR>::adjacencyFromVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    if(j < vertex) {
        const std::size_t vertexAdjacent = j;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertexAdjacent, vertex);
        return Adjacency<>(vertexAdjacent, edgeAdjacent);
    }
    else {
        const std::size_t vertexAdjacent = j + 1;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertex, vertexAdjacent);
        return Adjacency<>(vertexAdjacent, edgeAdjacent);
    }
}

/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline Adjacency<>
CompleteGraph<VISITOR>::adjacencyToVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    return adjacencyFromVertex(vertex, j);
}

/// Search for an edge (in constant time).
///
/// Indexing: findEdge(vertex0, vertex1)
///    - 0 1 2
///    0 - 3 4
///    1 3 - 5
///    2 4 5 -
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \return if an edge from vertex0 to vertex1 exists, pair.first is true
///     and pair.second is the index of such an edge. if no edge from vertex0
///     to vertex1 exists, pair.first is false and pair.second is undefined.
///
template<typename VISITOR>
inline std::pair<bool, std::size_t>
CompleteGraph<VISITOR>::findEdge(
    const std::size_t vertex0,
    const std::size_t vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());
    if(vertex0 == vertex1) {
        return std::pair<bool, std::size_t>(false, 0);
    }
    else if(vertex0 < vertex1) {
        return std::pair<bool, std::size_t>(true, edgeOfStrictlyIncreasingPairOfVertices(vertex0, vertex1));
    }
    else {
        return std::pair<bool, std::size_t>(true, edgeOfStrictlyIncreasingPairOfVertices(vertex1, vertex0));
    }
}

/// Indicate if multiple edges are enabled.
///
/// \return false
///
template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::multipleEdgesEnabled() const {
    return false;
}

// private
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::edgeOfStrictlyIncreasingPairOfVertices(
    const std::size_t vertex0,
    const std::size_t vertex1
) const {
    assert(vertex1 < numberOfVertices());
    assert(vertex0 < vertex1);
    return (numberOfVertices() - 1) * vertex0 - vertex0 * (vertex0 + 1) / 2 + vertex1 - 1;
}

// \cond SUPPRESS_DOXYGEN

// implementation of AdjacencyIterator

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator()
:   graph_(0),
    vertex_(0),
    adjacencyIndex_(0),
    adjacency_()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph
)
:   graph_(&graph),
    vertex_(0),
    adjacencyIndex_(0),
    adjacency_()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const std::size_t vertex
)
:   graph_(&graph),
    vertex_(vertex),
    adjacencyIndex_(0),
    adjacency_()
{
    assert(vertex < graph.numberOfVertices());
}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const std::size_t vertex,
    const std::size_t adjacencyIndex
)
:   graph_(&graph),
    vertex_(vertex),
    adjacencyIndex_(adjacencyIndex),
    adjacency_()
{
    assert(vertex < graph.numberOfVertices());
    assert(adjacencyIndex <= graph.numberOfVertices());
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator+=(
    const difference_type d
) {
    adjacencyIndex_ += d;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator-=(
    const difference_type d
) {
    adjacencyIndex_ -= d;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator++() {
    ++adjacencyIndex_;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator--() {
    --adjacencyIndex_;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator++(int) {
    AdjacencyIterator copy = *this;
    ++adjacencyIndex_;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator--(int) {
    AdjacencyIterator copy = *this;
    --adjacencyIndex_;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator+(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy += d;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator-(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy -= d;
    return copy;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator==(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ == other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator!=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ != other.adjacencyIndex_
        || vertex_ != other.vertex_
        || graph_ != other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator<(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ < other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator<=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ <= other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}


template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator>(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ > other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator>=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ >= other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline const Adjacency<>&
CompleteGraph<VISITOR>::AdjacencyIterator::operator*() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return adjacency_;
}

template<typename VISITOR>
inline const Adjacency<>*
CompleteGraph<VISITOR>::AdjacencyIterator::operator->() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return &adjacency_;
}

template<typename VISITOR>
inline const Adjacency<>&
CompleteGraph<VISITOR>::AdjacencyIterator::operator[](
    const std::size_t j
) {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_ + j);
    return adjacency_;
}

// implementation of VertexIterator

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator()
:   Base()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph
)
:   Base(graph)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const std::size_t vertex
)
:   Base(graph, vertex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const std::size_t vertex,
    const std::size_t adjacencyIndex
)
:   Base(graph, vertex, adjacencyIndex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const VertexIterator& it
)
:   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const AdjacencyIterator& it
)
:   Base(it)
{}

template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::VertexIterator::operator*() const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::VertexIterator::operator[](
    const std::size_t j
) const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

// implementation of EdgeIterator

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator()
:   Base()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph
)
:   Base(graph)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const std::size_t vertex
)
:   Base(graph, vertex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const std::size_t vertex,
    const std::size_t adjacencyIndex
)
:   Base(graph, vertex, adjacencyIndex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const EdgeIterator& it
)
:   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const AdjacencyIterator& it
)
:   Base(it)
{}

template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::EdgeIterator::operator*() const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::EdgeIterator::operator[](
    const std::size_t j
) const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

// \endcond

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_HXX
