#include <stdexcept>

#include "andres/graph/graph.hxx"
#include "andres/graph/bfs.hxx"

inline void test(const bool condition) {
    if(!condition) throw std::logic_error("test failed.");
}

struct SimpleCallback {
    SimpleCallback()
        :   vertexIndices_(),
            depths_()
        {}
    void operator()(
        const std::size_t v,
        const std::size_t depth,
        bool& proceed,
        bool& add
    ) {
        vertexIndices_.push_back(v);
        depths_.push_back(depth);
        proceed = true;
        add = true;
    }

    std::vector<std::size_t> vertexIndices_;
    std::vector<std::size_t> depths_;
};

struct SearchCallback {
    SearchCallback(const std::size_t vertexIndex)
        :   vertexIndex_(vertexIndex),
            vertexIndices_(),
            depths_()
        {}
    void operator()(
        const std::size_t v,
        const std::size_t depth,
        bool& proceed,
        bool& add
    ) {
        vertexIndices_.push_back(v);
        depths_.push_back(depth);
        proceed = (v != vertexIndex_);
        add = true;
    }

    std::size_t vertexIndex_;
    std::vector<std::size_t> vertexIndices_;
    std::vector<std::size_t> depths_;
};

struct BlockingCallback {
    BlockingCallback(const std::size_t vertexIndex)
        :   vertexIndex_(vertexIndex),
            vertexIndices_(),
            depths_()
        {}
    void operator()(
        const std::size_t v,
        const std::size_t depth,
        bool& proceed,
        bool& add
    ) {
        vertexIndices_.push_back(v);
        depths_.push_back(depth);
        proceed = true;
        add = (v != vertexIndex_);
    }

    std::size_t vertexIndex_;
    std::vector<std::size_t> vertexIndices_;
    std::vector<std::size_t> depths_;
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
        andres::graph::breadthFirstSearch(g, 0, callback);
        test(callback.vertexIndices_.size() == 7);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 1);
        test(callback.vertexIndices_[2] == 2);
        test(callback.vertexIndices_[3] == 4);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 5);
        test(callback.vertexIndices_[6] == 6);
        test(callback.depths_[0] == 0);
        test(callback.depths_[1] == 1);
        test(callback.depths_[2] == 1);
        test(callback.depths_[3] == 1);
        test(callback.depths_[4] == 2);
        test(callback.depths_[5] == 2);
        test(callback.depths_[6] == 2);
    }

    {
        SearchCallback callback(2);
        andres::graph::breadthFirstSearch(g, 0, callback);
        test(callback.vertexIndices_.size() == 3);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 1);
        test(callback.vertexIndices_[2] == 2);
        test(callback.depths_[0] == 0);
        test(callback.depths_[1] == 1);
        test(callback.depths_[2] == 1);
    }

    {
        BlockingCallback callback(2);
        andres::graph::breadthFirstSearch(g, 0, callback);
        test(callback.vertexIndices_.size() == 6);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 1);
        test(callback.vertexIndices_[2] == 2);
        test(callback.vertexIndices_[3] == 4);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 5);
        test(callback.depths_[0] == 0);
        test(callback.depths_[1] == 1);
        test(callback.depths_[2] == 1);
        test(callback.depths_[3] == 1);
        test(callback.depths_[4] == 2);
        test(callback.depths_[5] == 2);
    }

    // disconnected graph
    g.insertVertices(2);
    g.insertEdge(7, 8);

    {
        SimpleCallback callback;
        andres::graph::breadthFirstSearch(g, 0, callback);
        andres::graph::breadthFirstSearch(g, 7, callback);
        test(callback.vertexIndices_.size() == 9);
        test(callback.vertexIndices_[0] == 0);
        test(callback.vertexIndices_[1] == 1);
        test(callback.vertexIndices_[2] == 2);
        test(callback.vertexIndices_[3] == 4);
        test(callback.vertexIndices_[4] == 3);
        test(callback.vertexIndices_[5] == 5);
        test(callback.vertexIndices_[6] == 6);
        test(callback.vertexIndices_[7] == 7);
        test(callback.vertexIndices_[8] == 8);
        test(callback.depths_[0] == 0);
        test(callback.depths_[1] == 1);
        test(callback.depths_[2] == 1);
        test(callback.depths_[3] == 1);
        test(callback.depths_[4] == 2);
        test(callback.depths_[5] == 2);
        test(callback.depths_[6] == 2);
        test(callback.depths_[7] == 0);
        test(callback.depths_[8] == 1);
    }

    return 0;
}
