#pragma once
#ifndef ANDRES_GRAPH_FAST_MARCHING_HXX
#define ANDRES_GRAPH_FAST_MARCHING_HXX

#include <stdexcept>
#include <iterator>
#include <vector>
#include <type_traits>
        
#include "grid-graph.hxx"

namespace andres {
namespace graph {

template<class T = double, class S = std::size_t>
class FastMarchingBuffers {
public:
    typedef T value_type;
    typedef S size_type;

    FastMarchingBuffers(std::size_t N)
        :   targetEdgeOfVertexWithSource_(N),
            isVertexFrozen_(N),
            distances_(N)
        {}

    std::vector<size_type> targetEdgeOfVertexWithSource_;
    std::vector<unsigned char> isVertexFrozen_;
    std::vector<value_type> distances_;
};

template<typename T, typename Enable = void>
struct my_make_signed {
    typedef typename std::make_signed<T>::type type;
};

template<typename T>
struct my_make_signed<T,
    typename std::enable_if<std::is_floating_point<T>::value>::type> {
    typedef T type;
};

/// Fast Marching for (a subgraph of) a 2-dimensional grid graph.
///
/// Implementation by Margret Keuper <keuper@mpi-inf.mpg.de>
/// of the algorithm defined in:
///
/// J.A. Sethian. Level Set Methods and Fast Marching Methods. In: 
/// Evolving Interfaces in Computational Geometry, Fluid Mechanics, 
/// Computer Vision and Materials Science. Cambridge University Press, 
/// 1999.
///
/// \param inputGraph input graph.
/// \param inputEdgeValues iterator to the beginning of a sequence of edge values
/// \param sourceVertex source vertex.
/// \param targetGraph graph with the same vertices as inputGraph but, possibly, additional edges for which distances are to be computed.
/// \param interpolationOrder 0=nearest neighbor, 1=linear interpolation.
/// \param targetEdgeValues input AND output: input: infinity=to be computed, not infinity=known and fixed.
///
template<
    class INPUT_GRAPH_VISITOR,
    class INPUT_EDGE_VALUE_ITERATOR,
    class TARGET_GRAPH,
    class TARGET_EDGE_VALUE_ITERATOR
>
void
fastMarching(
    const andres::graph::GridGraph<2, INPUT_GRAPH_VISITOR>& inputGraph,
    INPUT_EDGE_VALUE_ITERATOR inputEdgeValues,
    const std::size_t sourceVertex,
    const TARGET_GRAPH& targetGraph,
    TARGET_EDGE_VALUE_ITERATOR targetEdgeValues,
    const std::size_t interpolationOrder
) {
    auto buffers = FastMarchingBuffers<typename TARGET_EDGE_VALUE_ITERATOR::value_type, std::size_t>();
    fastMarching(inputGraph, inputEdgeValues, sourceVertex, targetGraph, targetEdgeValues, interpolationOrder, buffers);
}

template<
    class INPUT_GRAPH_VISITOR,
    class INPUT_EDGE_VALUE_ITERATOR,
    class TARGET_GRAPH,
    class TARGET_EDGE_VALUE_ITERATOR
>
void
fastMarching(
    const andres::graph::GridGraph<2, INPUT_GRAPH_VISITOR>& inputGraph,
    INPUT_EDGE_VALUE_ITERATOR inputEdgeValues,
    const std::size_t sourceVertex,
    const TARGET_GRAPH& targetGraph,
    TARGET_EDGE_VALUE_ITERATOR targetEdgeValues,
    const std::size_t interpolationOrder,
    FastMarchingBuffers<typename TARGET_EDGE_VALUE_ITERATOR::value_type, std::size_t>& buffers
) {
    typedef std::size_t size_type;
    typedef andres::graph::GridGraph<2, INPUT_GRAPH_VISITOR> InputGraph;
    typedef typename InputGraph::AdjacencyType Adjacency;
    typedef typename InputGraph::EdgeCoordinate EdgeCoordinate;
    typedef typename std::iterator_traits<TARGET_EDGE_VALUE_ITERATOR>::value_type Value;
    typedef typename my_make_signed<Value>::type signedValue;
    if(interpolationOrder < 0 || interpolationOrder > 1) {
        throw std::runtime_error("specified interpolation order not implemented.");
    }

    if(inputGraph.numberOfVertices() != targetGraph.numberOfVertices()) {
        throw std::runtime_error("inputGraph.numberOfVertices() != targetGraph.numberOfVertices().");
    }

    const Value infinity = std::numeric_limits<Value>::has_infinity
        ? std::numeric_limits<Value>::infinity()
        : std::numeric_limits<Value>::max();
    
    // make explicit, for each vertex:
    // - if it is a neighbor of sourceVertex in the targret graph
    // - if so, what the index of the connecting edge is (-1 if none) 
    // initialize isVertexFrozen is a Boolean vector with
    // - size targetGraph.numberOfVertices()
    // - isVertexFrozen[v] == 1 iff, for the edge e = (sourceVertex, v) in the target graph, targetEdgeValues[e] < infinity
    size_type maxEdge = targetGraph.numberOfEdges();

    // buffers.targetEdgeOfVertexWithSource_.resize(targetGraph.numberOfVertices());
    std::fill(
        buffers.targetEdgeOfVertexWithSource_.begin(), 
        buffers.targetEdgeOfVertexWithSource_.end(), 
        maxEdge + 1
    );
    
    size_type numberOfVerticesLeft = targetGraph.numberOfEdgesFromVertex(sourceVertex);

    // buffers.isVertexFrozen_.resize(targetGraph.numberOfVertices());
    std::fill(buffers.isVertexFrozen_.begin(), buffers.isVertexFrozen_.end(), 0);
    
    // buffers.distances_.resize(targetGraph.numberOfVertices());
    std::fill(buffers.distances_.begin(), buffers.distances_.end(), infinity);
    
    buffers.distances_[sourceVertex] = 0;
 
    for (auto it = targetGraph.adjacenciesFromVertexBegin(sourceVertex); it != targetGraph.adjacenciesFromVertexEnd(sourceVertex); ++it)
        buffers.targetEdgeOfVertexWithSource_[it->vertex()] = it->edge();
    
    struct elem
    {
        size_type v;
        Value minArrivalTime;

        elem(size_type _v, Value _minArrivalTime) :
            v(_v), minArrivalTime(_minArrivalTime)
        {}

        // because we need min-priority queue
        bool operator<(elem const& other) const
        {
            return minArrivalTime > other.minArrivalTime;
        }
    };

    std::priority_queue<elem> trial;
    trial.emplace(sourceVertex, Value());

    while (!trial.empty() && numberOfVerticesLeft > 0)
    {
        while (!trial.empty() && buffers.isVertexFrozen_[trial.top().v])
            trial.pop();

        if (trial.empty())
            return;

        auto frozenVertex = trial.top().v;
        auto minArrivalTime = trial.top().minArrivalTime;
        
        trial.pop();

        buffers.isVertexFrozen_[frozenVertex] = 1;

        const size_type targetEdge = buffers.targetEdgeOfVertexWithSource_[frozenVertex];
        if (targetEdge <= maxEdge) // if frozenVertex is a neighbor of sourceVertex in targetGraph
        {
#pragma omp critical
            {
                targetEdgeValues[targetEdge] = minArrivalTime;
            }

            --numberOfVerticesLeft;
            
            if (numberOfVerticesLeft == 0)
                return; // all target edge values of interest computed (exactly)
        }
        
        // For each neighbor u of frozenVertex do unless isFrozen[u]
        // 1. Calculate arrival-time
        // 2. Insert u into trial set or update arrival time of u if already in the trial set.
        size_type numAdjacenciesOfFrozen = inputGraph.numberOfEdgesFromVertex(frozenVertex);
        for (size_type a_id = 0; a_id < numAdjacenciesOfFrozen; ++a_id)
        {
            const Adjacency& a = inputGraph.adjacencyFromVertex(frozenVertex, a_id);
            const size_type u = a.vertex();
            const size_type ue = a.edge();

            if (buffers.isVertexFrozen_[u] == 0)
            {
                // Compute Arrival time
                std::array<Value, 2> d;
                std::array<Value, 2> T = {{infinity,infinity}};
                
                EdgeCoordinate veCoord;
                for (size_type au_id = 0; au_id < inputGraph.numberOfEdgesFromVertex(u); ++au_id)
                {
                    const Adjacency& au = inputGraph.adjacencyFromVertex(u, au_id);
                    const size_type v = au.vertex();
                    const size_type ve = au.edge();
                    
                    inputGraph.edge(ve, veCoord);

                    const size_type dimension = veCoord.dimension_; // Dimension along which the edge runs (0=horiz)
                    const Value edgeWeight = inputEdgeValues[ve];
                    const Value vDistance = buffers.distances_[v];
                    
                    if (vDistance == infinity)
                        continue; // Otherwise addition overflows

                    if (T[dimension] == infinity || vDistance + edgeWeight < T[dimension] + d[dimension])
                    {
                        d[dimension] = edgeWeight;
                        T[dimension] = vDistance;
                    }
                }
                // Compute arrival time.


                Value at;
                if (T[0] < infinity)
                    if (T[1] < infinity)
                        at = std::min(T[1] + d[1], T[0] + d[0]);
                    else 
                        at = T[0] + d[0];
                else
                    if (T[1] < infinity) 
                        at = T[1] + d[1];
                    else
                        at = T[1];//set to infinity in this case
        
                if (interpolationOrder == 1)
                {
                    const Value inverseSpeed = 1;  // unexposed algorithm parameter. if large causality is violated.
                    if ((T[0] < infinity) && (T[1] < infinity) && (T[0] + d[0] > T[1]) && (T[1] + d[1] > T[0]))
                    {
                        const signedValue dx2 = d[0]*d[0];
                        const signedValue dy2 = d[1]*d[1];
                        const signedValue denom = dx2 + dy2;
                        const signedValue summand = (dy2*T[0] + dx2*T[1])/denom;
                        const signedValue nom = inverseSpeed*dx2*dy2 - dy2*T[0]*T[0] - dx2*T[1]*T[1];
                        const signedValue squared = nom/denom + summand*summand;
                        
                        if(squared >= 0)
                            at = std::min(at, static_cast<Value>(summand + std::sqrt(squared)));
                    }
                }

                buffers.distances_[u] = std::min(buffers.distances_[u], at);

                trial.emplace(u, buffers.distances_[u]);
            }
        }
    }
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_FAST_MARCHING_HXX
