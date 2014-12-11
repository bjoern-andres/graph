#pragma once
#ifndef ANDRES_GRAPH_BFS_HXX
#define ANDRES_GRAPH_BFS_HXX

#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>
#include <queue>

namespace andres {
namespace graph {

template<class S = std::size_t>
class BreadthFirstSearchData {
public:
    typedef S size_type;

    static const size_type NOT_VISITED;

    BreadthFirstSearchData(const size_type size = 0)
        :   depth_(size, NOT_VISITED),
            queue_()
        {}
    void assign(const size_type size = 0)
        {
            depth_.assign(size, NOT_VISITED);
            if(!queue_.empty()) {
                queue_ = std::queue<size_type>();
            }
        }
    size_type add(const size_type v, const size_type depth)
        { depth_[v] = depth; queue_.push(v); }
    size_type next()
        { const size_type v = queue_.front(); queue_.pop(); return v; }
    bool isEmpty() const
        { return queue_.empty(); }
    size_type depth(const size_type v) const
        { return depth_[v]; }

private:
    std::vector<size_type> depth_;
    std::queue<size_type> queue_;
};

template<class S>
const S BreadthFirstSearchData<S>::NOT_VISITED = std::numeric_limits<S>::max();

template<class GRAPH, class CALLBACK>
void
breadthFirstSearch(
    const GRAPH& g,
    const std::size_t startVertex,
    CALLBACK& callback,
    BreadthFirstSearchData<std::size_t>& data
) {
    typedef std::size_t size_type;
    typedef typename GRAPH::VertexIterator VertexIterator;
    typedef BreadthFirstSearchData<std::size_t> SearchDataType;

    assert(startVertex < g.numberOfVertices());
    data.assign(g.numberOfVertices());
    {
        bool proceed;
        bool add;
        callback(startVertex, 0, proceed, add);
        if(!proceed) {
            return;
        }
        if(add) {
            data.add(startVertex, 0);
        }
    }
    while(!data.isEmpty()) {
        const size_type v = data.next();
        const size_type depth = data.depth(v) + 1;
        for(VertexIterator it = g.verticesFromVertexBegin(v); it != g.verticesFromVertexEnd(v); ++it) {
            const size_type w = *it;
            if(data.depth(w) == SearchDataType::NOT_VISITED) {
                bool proceed;
                bool add;
                callback(w, depth, proceed, add);
                if(!proceed) {
                    return;
                }
                if(add) {
                    data.add(w, depth);
                }
            }
        }
    }
}

template<class GRAPH, class CALLBACK>
void
breadthFirstSearch(
    const GRAPH& g,
    CALLBACK& callback,
    BreadthFirstSearchData<std::size_t>& data
) {
    typedef BreadthFirstSearchData<std::size_t> SearchDataType;

    breadthFirstSearch(g, 0, callback, data);
    for(std::size_t v = 0; v < g.numberOfVertices(); ++v) {
        if(data.depth(v) == SearchDataType::NOT_VISITED) {
            breadthFirstSearch(g, v, callback, data);
        }
    }
}

template<class GRAPH, class CALLBACK>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const std::size_t startVertex,
    CALLBACK& callback
) {
    BreadthFirstSearchData<std::size_t> data;
    breadthFirstSearch(g, startVertex, callback, data);
}

template<class GRAPH, class CALLBACK>
inline void
breadthFirstSearch(
    const GRAPH& g,
    CALLBACK& callback
) {
    BreadthFirstSearchData<std::size_t> data;
    breadthFirstSearch(g, callback, data);
}

} // namespace graph
} // namespace andres

#endif
