#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/dfs.hxx"

inline void test(const bool condition) {
    if(!condition) throw std::logic_error("test failed.");
}

struct SimpleCallback {
    SimpleCallback()
        :   vertexIndices_()
        {}
    void operator()(
        const std::size_t v,
        bool& proceed,
        bool& addNeighbors
    ) {
        vertexIndices_.push_back(v);
        proceed = true;
        addNeighbors = true;
    }

    std::vector<std::size_t> vertexIndices_;
};

struct SearchCallback {
    SearchCallback(const std::size_t vertexIndex)
        :   vertexIndex_(vertexIndex),
            vertexIndices_()
        {}
    void operator()(
        const std::size_t v,
        bool& proceed,
        bool& addNeighbors
    ) {
        vertexIndices_.push_back(v);
        proceed = (v != vertexIndex_);
        addNeighbors = true;
    }

    std::size_t vertexIndex_;
    std::vector<std::size_t> vertexIndices_;
};

struct BlockingCallback {
    BlockingCallback(const std::size_t vertexIndex)
        :   vertexIndex_(vertexIndex),
            vertexIndices_()
        {}
    void operator()(
        const std::size_t v,
        bool& proceed,
        bool& addNeighbors
    ) {
        vertexIndices_.push_back(v);
        proceed = true;
        addNeighbors = (v != vertexIndex_);
    }

    std::size_t vertexIndex_;
    std::vector<std::size_t> vertexIndices_;
};

int main() {
    andres::graph::Graph<> g;
    g.insertVertices(7);
    g.insertEdge(0, 1);
    g.insertEdge(0, 2);
    g.insertEdge(0, 4);
    g.insertEdge(1, 3);
    g.insertEdge(1, 5);
    g.insertEdge(2, 6);
    g.insertEdge(4, 5);

    {
        SimpleCallback callback;
        andres::graph::depthFirstSearch(g, 0, callback);
        test(callback.vertexIndices_.size() == 7);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 4);
        test(callback.vertexIndices_[2] == 5);
        test(callback.vertexIndices_[3] == 1);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 2);
        test(callback.vertexIndices_[6] == 6);
    }

    {
        SearchCallback callback(2);
        andres::graph::depthFirstSearch(g, 0, callback);
        test(callback.vertexIndices_.size() == 6);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 4);
        test(callback.vertexIndices_[2] == 5);
        test(callback.vertexIndices_[3] == 1);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 2);
    }

    {
        BlockingCallback callback(2);
        andres::graph::depthFirstSearch(g, 0, callback);
        test(callback.vertexIndices_.size() == 6);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 4);
        test(callback.vertexIndices_[2] == 5);
        test(callback.vertexIndices_[3] == 1);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 2);
    }

    // disconnected graph
    g.insertVertices(2);
    g.insertEdge(7, 8);

    {
        SimpleCallback callback;
        andres::graph::depthFirstSearch(g, callback);
        test(callback.vertexIndices_.size() == 9);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 4);
        test(callback.vertexIndices_[2] == 5);
        test(callback.vertexIndices_[3] == 1);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 2);
        test(callback.vertexIndices_[6] == 6);
        test(callback.vertexIndices_[7] == 7);
        test(callback.vertexIndices_[8] == 8);
    }

    return 0;
}
