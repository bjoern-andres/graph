#pragma once
#ifndef UTILS_HXX
#define UTILS_HXX

#include <stack>
#include <vector>
#include <andres/graph/dfs.hxx>

template<typename GRAPH, typename VLA, typename ELA>
inline
void vertexToEdgeLabels(const GRAPH& graph, const VLA& vertex_labels, ELA& edge_labels)
{
    for (std::size_t i = 0; i < graph.numberOfEdges(); ++i)
    {
        auto v0 = graph.vertexOfEdge(i, 0);
        auto v1 = graph.vertexOfEdge(i, 1);

        edge_labels[i] = vertex_labels[v0] == vertex_labels[v1] ? 0 : 1;
    }
}

template<typename ORIGGRAPH, typename LIFTGRAPH, typename VLA, typename ELA>
inline
void vertexToEdgeLabels(const ORIGGRAPH& original_graph, const LIFTGRAPH& lifted_graph, const VLA& vertex_labels, ELA& edge_labels)
{
    std::fill(edge_labels.begin(), edge_labels.end(), 1);

    std::stack<typename VLA::value_type> S;
    std::vector<char> visited(original_graph.numberOfVertices());
    VLA new_labels(original_graph.numberOfVertices());

    // coonectedness on the original graph
    for (std::size_t i = 0, cnt = 0; i < original_graph.numberOfVertices(); ++i)
        if (!visited[i])
        {
            S.push(i);
            visited[i] = 1;

            auto label = vertex_labels[i];

            new_labels[i] = cnt;

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                    if (vertex_labels[it->vertex()] == label && !visited[it->vertex()])
                    {
                        S.push(it->vertex());
                        visited[it->vertex()] = 1;
                        new_labels[it->vertex()] = cnt;
                    }
            }

            ++cnt;
        }

    vertexToEdgeLabels(lifted_graph, new_labels, edge_labels);
}

template<typename GRAPH, typename VLA, typename ELA>
inline
void edgeToVertexLabels(const GRAPH& graph, const ELA& edge_labels, VLA& vertex_labels)
{
    struct mask
    {
        mask(const ELA& edge_labels) : edge_labels_(edge_labels)
        {}
        bool vertex(std::size_t i) const
            { return true; }
        bool edge(std::size_t i) const
            { return !edge_labels_[i]; }

        const ELA& edge_labels_;
    };

    andres::graph::DepthFirstSearchData<> dfs_data(graph.numberOfVertices());

    for (std::size_t i = 0, label = 0; i < graph.numberOfVertices(); ++i)
        if (!dfs_data.visited(i))
        {
            depthFirstSearch(
                graph,
                mask(edge_labels),
                i,
                [&](std::size_t v, bool& proceed, bool& add_neighbors)
                {
                    vertex_labels[v] = label;
                    proceed = true;
                    add_neighbors = true;
                },
                dfs_data);

            ++label;
        }
}

#endif // #ifndef UTILS_HXX
