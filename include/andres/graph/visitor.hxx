#pragma once
#ifndef ANDRES_GRAPH_VISITOR_HXX
#define ANDRES_GRAPH_VISITOR_HXX

#include <iostream>

namespace andres {
namespace graph {

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

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_VISITOR_HXX
