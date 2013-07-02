#pragma once
#ifndef ANDRES_GRAPH_MAX_FLOW_HXX
#define ANDRES_GRAPH_MAX_FLOW_HXX

#include <cassert>
#include <vector>
#include <queue>
#include <limits> // std::numeric_limit
#include <algorithm> // std::max, std::min

#include "andres/graph/shortest-paths.hxx"

namespace andres {
namespace graph {

/// Push-Relabel Algorithm for computing the maximum s-t-flow of a Digraph.
///
/// With FIFO vertex selection rule.
///
template<class GRAPH, class FLOW>
class MaxFlowPushRelabel {
public:
    typedef GRAPH GraphType;
    typedef FLOW Flow;
    typedef typename GraphType::VertexIterator VertexIterator;
    typedef typename GraphType::EdgeIterator EdgeIterator;
    typedef typename GraphType::AdjacencyIterator AdjacencyIterator;

    MaxFlowPushRelabel();
    void clear();
    Flow maxFlow() const;
    Flow flow(const size_t) const;
    size_t numberOfPushes() const;
    size_t numberOfRelabels() const;
    template<class EDGE_WEIGHT_ITERATOR>
        MaxFlowPushRelabel(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);
    template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
        MaxFlowPushRelabel(const GraphType&, const SUBGRAPH_MASK&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);
    template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
        Flow operator()(const GraphType&, const SUBGRAPH_MASK&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);

private:
    template<class EDGE_WEIGHT_ITERATOR>
        void push(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t);
    template<class EDGE_WEIGHT_ITERATOR>
        void pushBack(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t);
    template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
        void relabel(const GraphType&, const SUBGRAPH_MASK&, EDGE_WEIGHT_ITERATOR, const size_t);
    template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
        void discharge(const GraphType&, const SUBGRAPH_MASK&, EDGE_WEIGHT_ITERATOR, const size_t);
    void gapRelabel(const GraphType&, const size_t);

    std::vector<size_t> height_;
    std::vector<size_t> labelCount_;    
    std::vector<Flow> excess_;
    std::vector<Flow> flow_;
    std::vector<bool> active_;
    std::queue<size_t> queue_;
    size_t sourceVertexIndex_;
    size_t sinkVertexIndex_;
    size_t pushCount_;
    size_t relabelCount_;
};

/// Construct an instance of the push-relabel algorithm.
///
template<class GRAPH, class FLOW>
inline
MaxFlowPushRelabel<GRAPH, FLOW>::MaxFlowPushRelabel()
:   height_(),
    labelCount_(),
    excess_(), 
    flow_(),
    active_(),
    queue_(),
    sourceVertexIndex_(),
    sinkVertexIndex_(),
    pushCount_(), 
    relabelCount_()
{}

/// Construct an instance of the push-relabel algorithm.
///
/// \param graph A graph.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param sourceVertexIndex Index of the source vertex.
/// \param sinkVertexIndex Index of the sink vertex.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline
MaxFlowPushRelabel<GRAPH, FLOW>::MaxFlowPushRelabel(
    const GraphType& graph,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t sourceVertexIndex,
    const size_t sinkVertexIndex
)
:   height_(), 
    labelCount_(),
    excess_(), 
    flow_(),
    active_(),
    queue_(),
    sourceVertexIndex_(),
    sinkVertexIndex_(),
    pushCount_(), 
    relabelCount_()
{
    (*this)(graph, DefaultSubgraphMask<>(), edgeWeightIterator, sourceVertexIndex, sinkVertexIndex);
}

/// Construct an instance of the push-relabel algorithm.
///
/// \param graph A graph.
/// \param mask A subgraph mask.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param sourceVertexIndex Index of the source vertex.
/// \param sinkVertexIndex Index of the sink vertex.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
inline
MaxFlowPushRelabel<GRAPH, FLOW>::MaxFlowPushRelabel(
    const GraphType& graph,
    const SUBGRAPH_MASK& mask,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t sourceVertexIndex,
    const size_t sinkVertexIndex
)
:   height_(),
    labelCount_(),
    excess_(), 
    flow_(),
    active_(),
    queue_(),
    sourceVertexIndex_(),
    sinkVertexIndex_(),
    pushCount_(), 
    relabelCount_()
{
    (*this)(graph, mask, edgeWeightIterator, sourceVertexIndex, sinkVertexIndex);
}
    
/// Clear all members.
///
template<class GRAPH, class FLOW>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::clear() {
    height_.clear();
    labelCount_.clear();
    excess_.clear();
    flow_.clear();
    active_.clear();
    sourceVertexIndex_ = 0;
    sinkVertexIndex_ = 0;
    pushCount_ = 0;
    relabelCount_ = 0;

    assert(queue_.empty());
}

/// Return the maximum flow through the graph.
///
/// \return The excess flow at the sink vertex, which is equal to the max flow.
///
template<class GRAPH, class FLOW>
inline typename MaxFlowPushRelabel<GRAPH, FLOW>::Flow 
MaxFlowPushRelabel<GRAPH, FLOW>::maxFlow() const  {
    return excess_[sinkVertexIndex_];
}

/// Return the flow in a certain edge.
///
/// \param edgeIndex Index of an edge.
/// \return The flow in the edge given by edgeIndex.
///
template<class GRAPH, class FLOW>
inline typename MaxFlowPushRelabel<GRAPH, FLOW>::Flow 
MaxFlowPushRelabel<GRAPH, FLOW>::flow(
    const size_t edgeIndex
) const {
    assert(edgeIndex < flow_.size());
    return flow_[edgeIndex];
}

/// Return the total number of pushes executed.
///
/// \return The number of pushes.
///
template<class GRAPH, class FLOW>
inline size_t 
MaxFlowPushRelabel<GRAPH, FLOW>::numberOfPushes() const {
    return pushCount_;
}

/// Return the total number of relabels executed.
///
/// \return The number of relabels.
///
template<class GRAPH, class FLOW>
inline size_t 
MaxFlowPushRelabel<GRAPH, FLOW>::numberOfRelabels() const {
    return relabelCount_;
}

/// Initialize members and executes push-relabel algorithm.
///
/// \param graph A graph.
/// \param mask A subgraph mask.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param sourceVertexIndex Index of the source vertex.
/// \param sinkVertexIndex Index of the sink vertex.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
inline typename MaxFlowPushRelabel<GRAPH, FLOW>::Flow 
MaxFlowPushRelabel<GRAPH, FLOW>::operator()(
    const GraphType& graph,
    const SUBGRAPH_MASK& mask,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t sourceVertexIndex,
    const size_t sinkVertexIndex
) {
    assert(mask.vertex(sourceVertexIndex) && mask.vertex(sinkVertexIndex));
    
    const size_t numberOfVertices = graph.numberOfVertices();
    const size_t numberOfEdges = graph.numberOfEdges();

    // initialization
    sourceVertexIndex_ = sourceVertexIndex;
    sinkVertexIndex_ = sinkVertexIndex;
    pushCount_ = 0;
    relabelCount_ = 0;
    height_.resize(numberOfVertices);
    excess_.resize(numberOfVertices);
    active_.resize(numberOfVertices);
    std::fill(height_.begin(), height_.end(), size_t());
    std::fill(excess_.begin(), excess_.end(), Flow());
    std::fill(active_.begin(), active_.end(), false);
    height_[sourceVertexIndex] = numberOfVertices;
    excess_[sourceVertexIndex] = std::numeric_limits<Flow>::max(); // this is supposed to be infinite
    active_[sourceVertexIndex] = true;
    active_[sinkVertexIndex] = true;
    labelCount_.resize((2 * numberOfVertices) + 2); // 2n + 1 is the maximum possible height for a vertex
    std::fill(labelCount_.begin(), labelCount_.end(), size_t());
    labelCount_[0] = numberOfVertices - 1;
    labelCount_[numberOfVertices] = 1;
    flow_.resize(numberOfEdges);
    std::fill(flow_.begin(), flow_.end(), Flow());

    // first, push as much flow as possible from the source to all adjacent vertices
    for(EdgeIterator it = graph.edgesFromVertexBegin(sourceVertexIndex); it != graph.edgesFromVertexEnd(sourceVertexIndex); ++it) {
        const size_t edgeIndex = *it;
        if (mask.edge(edgeIndex)) {
            const size_t v = graph.vertexOfEdge(edgeIndex, 1);
            if (mask.vertex(v)) {
                push(graph, edgeWeightIterator, edgeIndex);
            }
        }
    }
    
    while(!queue_.empty()) { // main loop
        const size_t u = queue_.front();
        queue_.pop();
        active_[u] = false;
        discharge(graph, mask, edgeWeightIterator, u);
        active_[u] = false; // TODO: why does putting active_[u] = false twice decrease the number of pushes?
    }
    
    return maxFlow();
}

/// Push flow forward along an outgoing edge from a node to its child node.
///
/// \param graph A graph.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param edgeIndex Index of an edge.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::push(
    const GRAPH& graph,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t edgeIndex
) {
    const size_t u = graph.vertexOfEdge(edgeIndex, 0);
    const size_t v = graph.vertexOfEdge(edgeIndex, 1);
    // push from u to v
    Flow pushAmount = std::min(excess_[u], edgeWeightIterator[edgeIndex] - flow_[edgeIndex]);
    flow_[edgeIndex] += pushAmount;
    excess_[u] -= pushAmount;
    excess_[v] += pushAmount;
    if (!active_[v] && excess_[v] > 0) {
        active_[v] = true;
        queue_.push(v);
    }
    pushCount_++;
}

/// Push flow backward along an incoming edge from a node to its parent node.
///
/// \param graph A graph.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param edgeIndex Index of an edge.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::pushBack(
    const GRAPH& graph,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t edgeIndex
) {
    const size_t u = graph.vertexOfEdge(edgeIndex, 0);
    const size_t v = graph.vertexOfEdge(edgeIndex, 1);
    // push back from v to u
    Flow pushAmount = std::min(excess_[v], flow_[edgeIndex]);
    flow_[edgeIndex] -= pushAmount;
    excess_[u] += pushAmount;
    excess_[v] -= pushAmount;
    if (!active_[u] && excess_[u] > 0) {
        active_[u] = true;
        queue_.push(u);
    }
    pushCount_++;
}

/// Increase height of u to 1 greater than the minimum height of its neighbors.  Execute a gap relabel if a gap in heights exists.
///
/// \param graph A graph.
/// \param mask A subgraph mask.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param u Index of a vertex to relabel.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::relabel(
    const GRAPH& graph,
    const SUBGRAPH_MASK& mask,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t u
) {
    size_t minHeight = 2 * graph.numberOfVertices();
    const size_t oldHeight = height_[u];
    for(EdgeIterator it = graph.edgesFromVertexBegin(u); it != graph.edgesFromVertexEnd(u); ++it) {
        const size_t edgeIndex = *it;
        if (mask.edge(edgeIndex)) {
            const size_t v = graph.vertexOfEdge(edgeIndex, 1); // edge is (u, v)
            if (mask.vertex(v)) {
                if(edgeWeightIterator[edgeIndex] - flow_[edgeIndex] > 0) {
                    minHeight = std::min(minHeight, height_[v]);
                }
            }
        }
    }
    for(EdgeIterator it = graph.edgesToVertexBegin(u); it != graph.edgesToVertexEnd(u); ++it) {
        const size_t edgeIndex = *it;
        if (mask.edge(edgeIndex)) {
            const size_t v = graph.vertexOfEdge(edgeIndex, 0); // edge is (v, u)
            if (mask.vertex(v)) {
                if(flow_[edgeIndex] > 0) {
                    minHeight = std::min(minHeight, height_[v]);
                }
            }
        }
    }
    height_[u] = minHeight + 1;
    if (!active_[u] && excess_[u] > 0) {
        active_[u] = true;
        queue_.push(u);
    }
    relabelCount_++;

    // gap relabeling
    labelCount_[oldHeight]--;
    labelCount_[height_[u]]++;
    if (labelCount_[oldHeight] == 0) {
        gapRelabel(graph, oldHeight);
    }
}

/// While there is excess flow at u, try to push flow to its neighbors. If no push is available, relabel u.
///
/// \param graph A graph.
/// \param mask A subgraph mask.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param u Index of a vertex to discharge.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::discharge(
    const GRAPH& graph,
    const SUBGRAPH_MASK& mask,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t u
) {
    while(excess_[u] > 0) {
        for(EdgeIterator it = graph.edgesFromVertexBegin(u); it != graph.edgesFromVertexEnd(u); ++it) {
            const size_t edgeIndex = *it;
            if (mask.edge(edgeIndex)) {
                const size_t v = graph.vertexOfEdge(edgeIndex, 1); // edge is the pair (u, v)
                if (mask.vertex(v)) {
                    if(edgeWeightIterator[edgeIndex] - flow_[edgeIndex] > 0 && height_[u] > height_[v]) {
                        push(graph, edgeWeightIterator, edgeIndex);
                    }
                }
            }
        }
        for(EdgeIterator it = graph.edgesToVertexBegin(u); it != graph.edgesToVertexEnd(u); ++it) {
            const size_t edgeIndex = *it;
            if (mask.edge(edgeIndex)) {
                const size_t v = graph.vertexOfEdge(edgeIndex, 0); // edge is the pair (v, u)
                if (mask.vertex(v)) {
                    if(flow_[edgeIndex] > 0 && height_[u] > height_[v]) {
                        pushBack(graph, edgeWeightIterator, edgeIndex);
                    }
                }
            }
        }
        relabel(graph, mask, edgeWeightIterator, u);
    }
}

/// If there is a gap in heights, relabel all vertices with height above the gap to a height greater than the height of the source vertex.
///
/// \param graph A graph.
/// \param threshold The height at which a gap exists.
/// 
template<class GRAPH, class FLOW>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::gapRelabel(
    const GRAPH& graph,
    const size_t threshold
) {
    for(size_t i = 0; i < graph.numberOfVertices(); i++) {
        if(height_[i] > threshold) {
            height_[i] = std::max(height_[i], graph.numberOfVertices() + 1);;
        }
    }
}

/// Edmonds-Karp Algorithm for computing the maximum s-t-flow of a Digraph.
///
template<class GRAPH, class FLOW>
class MaxFlowEdmondsKarp {
public:
    typedef GRAPH GraphType;
    typedef FLOW Flow;
    
    MaxFlowEdmondsKarp();
    template <class EDGE_WEIGHT_ITERATOR>
        MaxFlowEdmondsKarp(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);
    template <class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
        MaxFlowEdmondsKarp(const GraphType&, const SUBGRAPH_MASK&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);
    void clear();
    Flow maxFlow() const;
    Flow flow(const size_t) const;
    template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
        Flow operator()(const GraphType&, const SUBGRAPH_MASK&, const EDGE_WEIGHT_ITERATOR, const size_t, const size_t);
    
private:
    std::deque<size_t> augmentingPath_;
    Digraph<> rgraph_;
    std::vector<Flow> flow_;
    size_t sourceVertexIndex_;
    size_t sinkVertexIndex_;
    Flow maxFlow_;
    
    template <class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
    class ResidualMask {
    public:
        ResidualMask(const size_t e, const std::vector<Flow>& flow, const EDGE_WEIGHT_ITERATOR& ew, const SUBGRAPH_MASK& sm)
            : edges(e), f(flow), edgeWeightIterator(ew), subgraphMask(sm)
            {}
        bool vertex (const size_t v) const
            { return subgraphMask.vertex(v); }
        bool edge (const size_t e) const
            { return capacity(e) > Flow() && inMask(e) ; }
        Flow capacity (const size_t e) const
            {
                if (e >= edges) {
                    return f[e - edges];
                } else {
                    return edgeWeightIterator[e] - f[e];
                }
            }

    private:
        const size_t edges;
        const std::vector<Flow>& f;
        const EDGE_WEIGHT_ITERATOR& edgeWeightIterator;
        const SUBGRAPH_MASK& subgraphMask;
        bool inMask (const size_t e) const
            {
                if(e >= edges) { return subgraphMask.edge(e - edges); }
                else { return subgraphMask.edge(e); }
            }
    };
};

/// Construct an instance of the Edmonds-Karp algorithm.
///
template<class GRAPH, class FLOW>
inline
MaxFlowEdmondsKarp<GRAPH, FLOW>::MaxFlowEdmondsKarp()
:   augmentingPath_(),
    rgraph_(),
    flow_(),
    sourceVertexIndex_(),
    sinkVertexIndex_(),
    maxFlow_()
{}

/// Construct an instance of the Edmonds-Karp algorithm.
///
/// \param graph A graph.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param sourceVertexIndex Index of the source vertex.
/// \param sinkVertexIndex Index of the sink vertex.
///
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline
MaxFlowEdmondsKarp<GRAPH, FLOW>::MaxFlowEdmondsKarp(
    const GraphType& graph,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t sourceVertexIndex,
    const size_t sinkVertexIndex
)
:   augmentingPath_(),
    rgraph_(),
    flow_(),
    sourceVertexIndex_(),
    sinkVertexIndex_(),
    maxFlow_()
{
   (*this)(graph, DefaultSubgraphMask<>(), edgeWeightIterator, sourceVertexIndex, sinkVertexIndex);
}
    
/// Construct an instance of the Edmonds-Karp algorithm.
///
/// \param graph A graph.
/// \param subgraphMask A subgraph mask.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param sourceVertexIndex Index of the source vertex.
/// \param sinkVertexIndex Index of the sink vertex.
///
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
inline
MaxFlowEdmondsKarp<GRAPH, FLOW>::MaxFlowEdmondsKarp(
    const GraphType& graph,
    const SUBGRAPH_MASK& subgraphMask,
    EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t sourceVertexIndex,
    const size_t sinkVertexIndex
)
:   augmentingPath_(),
    rgraph_(),
    flow_(),
    sourceVertexIndex_(),
    sinkVertexIndex_(),
    maxFlow_()
{
   (*this)(graph, subgraphMask, edgeWeightIterator, sourceVertexIndex, sinkVertexIndex);
}

/// Clear all members.
///
template<class GRAPH, class FLOW>
inline void
MaxFlowEdmondsKarp<GRAPH, FLOW>::clear() {
    augmentingPath_.clear();
    rgraph_.assign();
    flow_.clear();
    sourceVertexIndex_ = 0;
    sinkVertexIndex_ = 0;
    maxFlow_ = Flow();
}

/// Return the maximum flow through the graph.
///
/// \return The max flow.
///
template<class GRAPH, class FLOW>
inline typename MaxFlowEdmondsKarp<GRAPH, FLOW>::Flow
MaxFlowEdmondsKarp<GRAPH, FLOW>::maxFlow() const  {
    return maxFlow_;
}

/// Return the flow in a certain edge.
///
/// \param edgeIndex Index of an edge.
/// \return The flow in the edge given by edgeIndex.
///
template<class GRAPH, class FLOW>
inline typename MaxFlowEdmondsKarp<GRAPH, FLOW>::Flow
MaxFlowEdmondsKarp<GRAPH, FLOW>::flow(
    const size_t edgeIndex
) const {
    assert(edgeIndex < flow_.size());
    return flow_[edgeIndex];
}

/// Initialize members and executes Edmonds-Karp algorithm.
///
/// \param graph A graph.
/// \param subgraphMask A subgraph mask.
/// \param edgeWeightIterator Iterator to the beginning of a sequence of edge weights.
/// \param sourceVertexIndex Index of the source vertex.
/// \param sinkVertexIndex Index of the sink vertex.
/// 
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR, class SUBGRAPH_MASK>
typename MaxFlowEdmondsKarp<GRAPH, FLOW>::Flow
MaxFlowEdmondsKarp<GRAPH, FLOW>::operator()(
    const GraphType& graph,
    const SUBGRAPH_MASK& subgraphMask,
    const EDGE_WEIGHT_ITERATOR edgeWeightIterator,
    const size_t sourceVertexIndex,
    const size_t sinkVertexIndex
) {
    const size_t numberOfEdges = graph.numberOfEdges();
    
    // initialization
    sourceVertexIndex_ = sourceVertexIndex;
    sinkVertexIndex_ = sinkVertexIndex;
    flow_.resize(numberOfEdges);
    std::fill(flow_.begin(), flow_.end(), Flow());
    augmentingPath_.clear();
    rgraph_.assign(graph.numberOfVertices());
    for(size_t edge = 0; edge < numberOfEdges; ++edge) {
        rgraph_.insertEdge(graph.vertexOfEdge(edge, 0), graph.vertexOfEdge(edge, 1));
    }
    for(size_t edge = 0; edge < numberOfEdges; ++edge) {
        rgraph_.insertEdge(graph.vertexOfEdge(edge, 1), graph.vertexOfEdge(edge, 0));
    }
    ResidualMask<EDGE_WEIGHT_ITERATOR, SUBGRAPH_MASK> rm(numberOfEdges, flow_, edgeWeightIterator, subgraphMask);
    maxFlow_ = Flow();
    
    // while there's an augmenting path, augment flow along it. choose the shortest augmenting path by number of edges
    Flow minResidualCapacity;
    while(spspEdges(rgraph_, rm, sourceVertexIndex_, sinkVertexIndex_, augmentingPath_)) {
        minResidualCapacity = rm.capacity(augmentingPath_[0]);
        for(std::deque<size_t>::iterator it = augmentingPath_.begin(); it < augmentingPath_.end(); ++it) {
            if(rm.capacity(*it) < minResidualCapacity) minResidualCapacity = rm.capacity(*it);
        }
        for(std::deque<size_t>::iterator it = augmentingPath_.begin(); it < augmentingPath_.end(); ++it) {
            if(*it < numberOfEdges) flow_[*it] += minResidualCapacity;
            else flow_[*it - numberOfEdges] -= minResidualCapacity;
        }
    }
    
    // sum flow out of source to get the max flow
    for(size_t edge = 0; edge < numberOfEdges; ++edge) {
        if(graph.vertexOfEdge(edge, 0) == sourceVertexIndex_) {
            maxFlow_ += flow_[edge];
        }
    }
    
    return maxFlow_;
}
    
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MAX_FLOW_HXX
