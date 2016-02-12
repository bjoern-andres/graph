#pragma once
#ifndef PROBABILISTIC_LIFTING
#define PROBABILISTIC_LIFTING

#include <andres/graph/shortest-paths.hxx>

#include "fast-marching.hxx"

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

        andres::graph::sssp(
            sourceGraph,
            andres::graph::DefaultSubgraphMask<>(),
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
    const andres::graph::GridGraph<2, SOURCE_GRAPH_VISITOR>& sourceGraph,
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

#endif