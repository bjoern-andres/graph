#pragma once
#ifndef ANDRES_GRAPH_LIFTING_HXX
#define ANDRES_GRAPH_LIFTING_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <iterator> // std::iterator_traits
#include <algorithm> // std::fill
#include <vector>

#include "grid-graph.hxx"
#include "bfs.hxx"
#include "shortest-paths.hxx"
#include "fast-marching.hxx"

namespace andres {
namespace graph {

enum class LiftingMetric { PathLength, L2 };

/// Lift a graph.
template<class INPUT_GRAPH, class OUTPUT_GRAPH>
inline void
lift(
    const INPUT_GRAPH& inputGraph,
    OUTPUT_GRAPH& outputGraph,
    const std::size_t distanceUpperBound,
    const std::size_t distanceLowerBound = 0
) {
    typedef std::size_t size_type;

    if(outputGraph.numberOfVertices() != 0)
        throw std::runtime_error("output graph is not empty.");
    
    outputGraph.insertVertices(inputGraph.numberOfVertices());

    BreadthFirstSearchData<size_type> breadthFirstSearchData(inputGraph.numberOfVertices());
    std::vector<std::size_t> visited;
    std::vector<std::size_t> vertices;

    for (size_type v = 0; v < inputGraph.numberOfVertices(); ++v)
    {
        breadthFirstSearch(
            inputGraph,
            v,
            [&](size_type w, size_type depth, bool& proceed, bool& add)
            {
                proceed = true;
                add = false;
                
                if (depth <= distanceUpperBound)
                {
                    if (depth + 1 <= distanceUpperBound)
                    {
                        add = true;
                        visited.push_back(w);        
                    }

                    if (depth > distanceLowerBound)
                        vertices.push_back(w);
                }
            },
            breadthFirstSearchData
            );

        std::sort(vertices.begin(), vertices.end());

        for (auto w : vertices)
            outputGraph.insertEdge(v, w);

        vertices.clear();

        breadthFirstSearchData.depth(v) = BreadthFirstSearchData<size_type>::NOT_VISITED;
        for (auto w : visited)
            breadthFirstSearchData.depth(w) = BreadthFirstSearchData<size_type>::NOT_VISITED;
        
        visited.clear();
    }
}

/// Lift a grid graph - a faster implementation using the grid structure.
template<class INPUT_GRAPH_VISITOR, class OUTPUT_GRAPH>
inline void
lift(
    const GridGraph<2, INPUT_GRAPH_VISITOR>& inputGraph,
    OUTPUT_GRAPH& outputGraph,
    const std::size_t distanceUpperBound,
    const std::size_t distanceLowerBound = 0,
    LiftingMetric metric = LiftingMetric::PathLength
) {
    typedef GridGraph<2, INPUT_GRAPH_VISITOR> INPUT_GRAPH;

    typedef std::size_t size_type;
    typedef typename INPUT_GRAPH::VertexCoordinate VertexCoordinate;

    const size_type distanceUpperBoundSquared = distanceUpperBound * distanceUpperBound;
    const size_type distanceLowerBoundSquared = distanceLowerBound * distanceLowerBound;

    if(outputGraph.numberOfVertices() != 0)
        throw std::runtime_error("output graph is not empty.");

    outputGraph.insertVertices(inputGraph.numberOfVertices());

    VertexCoordinate cv;
    for(size_type v = 0; v < inputGraph.numberOfVertices(); ++v)
    {
        inputGraph.vertex(v, cv);

        // fill above portion of the window
        if (cv[1] > 0)
        {
            const std::size_t row0 = cv[1] < distanceUpperBound ? 0 : cv[1] - distanceUpperBound;
            std::size_t offsetY = 1;

            // We use yPlus = y+1 to avoid the 0-1 case. In this block y = yPlus-1.
            for(std::size_t yPlus = cv[1]; yPlus > row0; --yPlus, ++offsetY)
            {
                const std::size_t offsetX = (metric == LiftingMetric::PathLength) ?
                                            distanceUpperBound - offsetY:
                                            ::floor(::sqrt(distanceUpperBoundSquared - offsetY * offsetY));

                const std::size_t col0 = cv[0] < offsetX ? 0 : cv[0] - offsetX;
                std::size_t colN = cv[0] + offsetX;
                if(colN > inputGraph.shape(0) - 1)
                    colN = inputGraph.shape(0) - 1;
                
                for (std::size_t x = col0; x <= colN; ++x)
                {
                    if (metric == LiftingMetric::PathLength)
                    {
                        const std::size_t distance = ::abs(x - cv[0]) + ::abs(yPlus - 1 - cv[1]);

                        if (distance > distanceLowerBound)
                        {
                            const size_type w = inputGraph.vertex({{x, yPlus - 1}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                    else
                    {
                        const std::size_t sqaredDistance = (x - cv[0]) * (x - cv[0]) + (yPlus - 1 - cv[1]) * (yPlus - 1 - cv[1]);

                        if (sqaredDistance > distanceLowerBoundSquared)
                        {
                            const size_type w = inputGraph.vertex({{x, yPlus - 1}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                }
            }
        }

        // Add middle horizontal line, except for the central point
        {
            const std::size_t offsetX = distanceUpperBound;
            const std::size_t y = cv[1];
            const std::size_t col0 = cv[0] < offsetX ? 0 : cv[0] - offsetX;
            std::size_t colN = cv[0] + offsetX;

            if (colN > inputGraph.shape(0) - 1)
                colN = inputGraph.shape(0) - 1;
            
            if (cv[0] > distanceLowerBound)
                for (std::size_t x = col0; x <= cv[0] - distanceLowerBound - 1; ++x)
                {
                    const size_type& w = inputGraph.vertex({{x, y}});
                    outputGraph.insertEdge(v, w);
                }

            for (std::size_t x = cv[0] + distanceLowerBound + 1; x <= colN; ++x)
            {
                const size_type& w = inputGraph.vertex({{x, y}});
                outputGraph.insertEdge(v, w);
            }
        }

        // fill below window
        if (cv[1] < inputGraph.shape(1) - 1)
        {
            const std::size_t row0 = cv[1] + 1;
            std::size_t rowN = cv[1] + distanceUpperBound;
            
            if (cv[1] + distanceUpperBound > inputGraph.shape(1) - 1)
                rowN = inputGraph.shape(1) - 1;

            std::size_t offsetY = 1;
            for (std::size_t y = row0; y <= rowN; ++y, ++offsetY)
            {
                const std::size_t offsetX = (metric == LiftingMetric::PathLength) ?
                                            distanceUpperBound - offsetY :
                                            ::floor(::sqrt(distanceUpperBoundSquared - offsetY * offsetY));
                
                const std::size_t col0 = cv[0] < offsetX ? 0 : cv[0] - offsetX;
                std::size_t colN = cv[0] + offsetX;

                if (colN > inputGraph.shape(0) - 1)
                    colN = inputGraph.shape(0) - 1;

                for (std::size_t x = col0; x <= colN; ++x)
                {
                    if (metric == LiftingMetric::PathLength)
                    {
                        const std::size_t distance = ::abs(x - cv[0]) + ::abs(y - cv[1]);
                        
                        if (distance > distanceLowerBound)
                        {
                            const size_type w = inputGraph.vertex({{x, y}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                    else
                    {
                        const std::size_t sqaredDistance = (x - cv[0]) * (x - cv[0]) + (y - cv[1]) * (y - cv[1]);

                        if (sqaredDistance > distanceLowerBoundSquared)
                        {
                            const size_type w = inputGraph.vertex({{x, y}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                }
            }
        }
    }
} 
/// Lift edge values from a source graph to a target graph.
template<class SOURCE_GRAPH, class TARGET_GRAPH, class SOURCE_EDGE_VALUE_ITERATOR, class TARGET_EDGE_VALUE_ITERATOR>
void
liftEdgeValues(
    const SOURCE_GRAPH& sourceGraph,
    const TARGET_GRAPH& targetGraph,
    const SOURCE_EDGE_VALUE_ITERATOR sit,
    TARGET_EDGE_VALUE_ITERATOR tit
) {
    typedef std::size_t size_type;
    typedef typename std::iterator_traits<SOURCE_EDGE_VALUE_ITERATOR>::value_type value_type;

    if(sourceGraph.numberOfVertices() != targetGraph.numberOfVertices())
        throw std::runtime_error("number of vertices mismatch between source and target graph.");

    class SSSPVisitor
    {
    public:
        SSSPVisitor(std::vector<char>& target_vertices, std::size_t number_of_target_vertices) :
            target_vertices_(target_vertices), vertices_left_(number_of_target_vertices)
        { }

        bool operator()(typename std::vector<value_type>::iterator distances, typename std::vector<size_type>::iterator parents, std::size_t v)
        {
            if (target_vertices_[v])
            {
                target_vertices_[v] = 0;
                --vertices_left_;
            }

            return vertices_left_;
        }

    private:
        std::vector<char>& target_vertices_;
        std::size_t vertices_left_ { 0 };
    };

    std::vector<value_type> distances(sourceGraph.numberOfVertices());
    std::vector<size_type> parents(sourceGraph.numberOfVertices());
    std::vector<char> target_vertices(sourceGraph.numberOfVertices());
    std::vector<std::size_t> target_vertices_edges(sourceGraph.numberOfVertices());

#pragma omp parallel for firstprivate(distances, parents, target_vertices, target_vertices_edges)
    for(long int i = 0; i < targetGraph.numberOfVertices(); ++i)
    {
        std::fill(target_vertices.begin(), target_vertices.end(), char());

        for (auto it = targetGraph.adjacenciesFromVertexBegin(i); it != targetGraph.adjacenciesFromVertexEnd(i); ++it)
        {
            target_vertices[it->vertex()] = 1;
            target_vertices_edges[it->vertex()] = it->edge();
        }

        SSSPVisitor visitor(target_vertices, targetGraph.numberOfEdgesFromVertex(i));

        sssp(
            sourceGraph,
            DefaultSubgraphMask<>(),
            i,
            sit,
            distances.begin(), // buffer
            parents.begin(), // buffer
            visitor
        );

#pragma omp critical
        {
            for (auto it = targetGraph.adjacenciesFromVertexBegin(i); it != targetGraph.adjacenciesFromVertexEnd(i); ++it)
                tit[it->edge()] = distances[it->vertex()];
        }
    }

    // restore original edge weights
    for (std::size_t e = 0; e < sourceGraph.numberOfEdges(); ++e)
    {
        auto v0 = sourceGraph.vertexOfEdge(e, 0);
        auto v1 = sourceGraph.vertexOfEdge(e, 1);

        tit[targetGraph.findEdge(v0, v1).second] = sit[e];
    }
}

/// Lift edge values from a 2-dimensional grid graph to a target graph.
template<class SOURCE_GRAPH_VISITOR, class TARGET_GRAPH, class SOURCE_EDGE_VALUE_ITERATOR, class TARGET_EDGE_VALUE_ITERATOR>
void
liftEdgeValues(
    const GridGraph<2, SOURCE_GRAPH_VISITOR>& sourceGraph,
    const TARGET_GRAPH& targetGraph,
    SOURCE_EDGE_VALUE_ITERATOR sit,
    TARGET_EDGE_VALUE_ITERATOR tit,
    const std::size_t interpolationOrder = 0
) {
    typedef std::size_t size_type;
    typedef typename std::iterator_traits<TARGET_EDGE_VALUE_ITERATOR>::value_type target_value_type;

    if(sourceGraph.numberOfVertices() != targetGraph.numberOfVertices()) {
        throw std::runtime_error("number of vertices mismatch between source and target graph.");
    }

    std::size_t cnt = 0;
    FastMarchingBuffers<target_value_type, size_type> buffers(sourceGraph.numberOfVertices());

#pragma omp parallel for firstprivate(buffers)
    for(size_type v = 0; v < targetGraph.numberOfVertices(); ++v)
        fastMarching(
            sourceGraph,
            sit,
            v,
            targetGraph,
            tit,
            interpolationOrder,
            buffers
        );

    // restore original edge weights
    for (std::size_t e = 0; e < sourceGraph.numberOfEdges(); ++e)
    {
        auto v0 = sourceGraph.vertexOfEdge(e, 0);
        auto v1 = sourceGraph.vertexOfEdge(e, 1);

        tit[targetGraph.findEdge(v0, v1).second] = sit[e];
    }
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_LIFTING_HXX
