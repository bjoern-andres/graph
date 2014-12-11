#pragma once
#ifndef ANDRES_GRAPH_DFS_HXX
#define ANDRES_GRAPH_DFS_HXX

#include <cassert>
#include <cstddef>
#include <vector>
#include <stack>

namespace andres {
namespace graph {

template<class S = std::size_t>
class DepthFirstSearchData {
public:
    typedef S size_type;

    DepthFirstSearchData(const size_type size = 0)
        :   visited_(size, 0),
            stack_()
        {}
    void assign(const size_type size = 0)
        {
            visited_.assign(size, 0);
            if(!stack_.empty()) {
                stack_ = std::stack<size_type>();
            }
        }
    size_type add(const size_type v)
        { stack_.push(v); }
    size_type next()
        { const size_type v = stack_.top(); stack_.pop(); return v; }
    bool isEmpty() const
        { return stack_.empty(); }
    unsigned char& visited(const size_type v)
        { return visited_[v]; }
    unsigned char visited(const size_type v) const
        { return visited_[v]; }

private:
    std::vector<unsigned char> visited_;
    std::stack<size_type> stack_;
};

template<class GRAPH, class CALLBACK>
void
depthFirstSearch(
    const GRAPH& g,
    const std::size_t startVertex,
    CALLBACK& callback,
    DepthFirstSearchData<std::size_t>& data
) {
    typedef std::size_t size_type;
    typedef typename GRAPH::VertexIterator VertexIterator;
    typedef DepthFirstSearchData<std::size_t> SearchDataType;

    assert(startVertex < g.numberOfVertices());
    data.assign(g.numberOfVertices());
    data.add(startVertex);
    while(!data.isEmpty()) {
        const size_type v = data.next();
        if(data.visited(v) == 0) {
            data.visited(v) = 1;
            bool proceed;
            bool addNeighbors;
            callback(v, proceed, addNeighbors);
            if(!proceed) {
                return;
            }
            if(addNeighbors) {
                for(VertexIterator it = g.verticesFromVertexBegin(v); it != g.verticesFromVertexEnd(v); ++it) {
                    const size_type w = *it;
                    data.add(w);
                }
            }
        }
    }
}

template<class GRAPH, class CALLBACK>
void
depthFirstSearch(
    const GRAPH& g,
    CALLBACK& callback,
    DepthFirstSearchData<std::size_t>& data
) {
    typedef DepthFirstSearchData<std::size_t> SearchDataType;

    depthFirstSearch(g, 0, callback, data);
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        if(data.visited(v) == 0) {
            depthFirstSearch(g, v, callback, data);
        }
    }
}

template<class GRAPH, class CALLBACK>
inline void
depthFirstSearch(
    const GRAPH& g,
    const std::size_t startVertex,
    CALLBACK& callback
) {
    DepthFirstSearchData<std::size_t> data;
    depthFirstSearch(g, startVertex, callback, data);
}

template<class GRAPH, class CALLBACK>
inline void
depthFirstSearch(
    const GRAPH& g,
    CALLBACK& callback
) {
    DepthFirstSearchData<std::size_t> data;
    depthFirstSearch(g, callback, data);
}

} // namespace graph
} // namespace andres

#endif
