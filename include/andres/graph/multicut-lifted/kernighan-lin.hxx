#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_KERNIGHAN_LIN_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_KERNIGHAN_LIN_HXX

#include <iomanip>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <stack>

#include "andres/graph/components.hxx"
#include "andres/graph/twocut-lifted/kernighan-lin.hxx"

namespace andres {
namespace graph {
namespace multicut_lifted {

struct KernighanLinSettings
{
    std::size_t numberOfInnerIterations { std::numeric_limits<std::size_t>::max() };
    std::size_t numberOfOuterIterations { 100 };
    double epsilon { 1e-7 };
    bool verbose { true };
};

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA>
void kernighanLin(
    const ORIGINAL_GRAPH& original_graph,
    const LIFTED_GRAPH& lifted_graph,
    const ECA& edgeCosts,
    const ELA& inputLabels,
    ELA& outputLabels,
    const KernighanLinSettings settings = KernighanLinSettings());

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA, typename VIS>
void kernighanLin(
    const ORIGINAL_GRAPH& original_graph,
    const LIFTED_GRAPH& lifted_graph,
    const ECA& edgeCosts,
    const ELA& inputLabels,
    ELA& outputLabels,
    VIS& visitor,
    const KernighanLinSettings settings = KernighanLinSettings());





template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA>
void kernighanLin(
    const ORIGINAL_GRAPH& original_graph,
    const LIFTED_GRAPH& lifted_graph,
    const ECA& edgeCosts,
    const ELA& inputLabels,
    ELA& outputLabels,
    const KernighanLinSettings settings = KernighanLinSettings())
{
    struct Visitor
    {
        constexpr bool operator()(std::vector<std::size_t>& vertex_labels)
        {
            return true;
        }

        constexpr bool time_limit_exceeded()
        {
            return false;
        }
    } visitor;

    kernighanLin(original_graph, lifted_graph, edgeCosts, inputLabels, outputLabels, visitor, settings);
}

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA, typename VIS>
inline
void kernighanLin(
    const ORIGINAL_GRAPH& original_graph,
    const LIFTED_GRAPH& lifted_graph,
    const ECA& edgeCosts,
    const ELA& inputLabels,
    ELA& outputLabels,
    VIS& visitor,
    const KernighanLinSettings settings)
{
    struct SubgraphWithCut { // a subgraph with cut mask
        SubgraphWithCut(const ELA& labels, std::vector<std::size_t> const& edge_in_lifted_graph)
            : labels_(labels), edge_in_lifted_graph_(edge_in_lifted_graph)
            {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return labels_[edge_in_lifted_graph_[e]] == 0; }

        std::vector<std::size_t> const& edge_in_lifted_graph_;
        const ELA& labels_;
    };

    std::vector<std::size_t> edge_in_lifted_graph(original_graph.numberOfEdges());
    for (std::size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto v0 = original_graph.vertexOfEdge(i, 0);
        auto v1 = original_graph.vertexOfEdge(i, 1);

        edge_in_lifted_graph[i] = lifted_graph.findEdge(v0, v1).second;
    }

    // build decomposition based on the current multicut
    ComponentsBySearch<ORIGINAL_GRAPH> components;
    components.build(original_graph, SubgraphWithCut(inputLabels, edge_in_lifted_graph));

    double current_energy_value = .0;

    // find out how many connected components there are
    // check if the input multicut labeling is valid
    std::size_t numberOfComponents = 0;
    for(std::size_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
    {
        outputLabels[edge] = inputLabels[edge];

        auto v0 = lifted_graph.vertexOfEdge(edge, 0);
        auto v1 = lifted_graph.vertexOfEdge(edge, 1);

        numberOfComponents = std::max(numberOfComponents, std::max(components.labels_[v0], components.labels_[v1]));

        if (inputLabels[edge])
            current_energy_value += edgeCosts[edge];

        if (static_cast<bool>(inputLabels[edge]) != !components.areConnected(v0, v1))
            throw std::runtime_error("the input multicut labeling is invalid.");
    }

    ++numberOfComponents;

    std::vector<std::vector<std::size_t>> partitions(numberOfComponents);

    twocut_lifted::TwoCutSettings twocut_settings;
    twocut_settings.numberOfIterations = settings.numberOfInnerIterations;
    twocut_settings.epsilon = settings.epsilon;

    twocut_lifted::TwoCutBuffers<ORIGINAL_GRAPH> twocut_buffers(original_graph);

    // build partitions
    for (std::size_t i = 0; i < components.labels_.size(); ++i)
    {
        partitions[components.labels_[i]].push_back(i);
        twocut_buffers.vertex_labels[i] = components.labels_[i];
    }

    twocut_buffers.max_not_used_label = partitions.size();

    if (settings.verbose)
    {
        std::cout << "Starting energy: " << current_energy_value << std::endl;
        std::cout << std::setw(4) << "Iter" << std::setw(16) << "True decrease" << std::setw(15) << "Total decrease" << std::setw(15) << "Pair updates" << std::setw(15) << "New sets" << std::setw(15) << "Num. of sets\n";
    }

    auto last_good_vertex_labels = twocut_buffers.vertex_labels;

    std::vector<char> visited(original_graph.numberOfVertices());

    // 1 if i-th partitioned changed since last iteration, 0 otherwise
    std::vector<char> changed(numberOfComponents, 1);

    // interatively update bipartition in order to minimize the total cost of the multicut
    for (std::size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {        
        auto energy_decrease = .0;

        std::vector<std::unordered_set<std::size_t>> edges(numberOfComponents);
        for (std::size_t e = 0; e < original_graph.numberOfEdges(); ++e)
        {
            auto v0 = twocut_buffers.vertex_labels[original_graph.vertexOfEdge(e, 0)];
            auto v1 = twocut_buffers.vertex_labels[original_graph.vertexOfEdge(e, 1)];

            if (v0 != v1)
                edges[std::min(v0, v1)].insert(std::max(v0, v1));
        }

        for (std::size_t i = 0; i < numberOfComponents && !visitor.time_limit_exceeded(); ++i)
            if (!partitions[i].empty())
                for (auto j = edges[i].begin(); j != edges[i].end() && !visitor.time_limit_exceeded(); ++j)
                    if (!partitions[*j].empty() && (changed[*j] || changed[i]))
                    {
                        auto ret = twocut_lifted::kernighanLin(original_graph, lifted_graph, edgeCosts, partitions[i], partitions[*j], twocut_buffers, twocut_settings);

                        if (ret > settings.epsilon)
                            changed[i] = changed[*j] = 1;

                        energy_decrease += ret;

                        if (partitions[i].size() == 0)
                            break;
                    }

        auto ee = energy_decrease;

        // remove partitions that became empty after the previous step
        auto new_end = std::partition(partitions.begin(), partitions.end(), [](const std::vector<std::size_t>& s) { return !s.empty(); });
        partitions.resize(new_end - partitions.begin());

        // try to intoduce new partitions
        for (std::size_t i = 0, p_size = partitions.size(); i < p_size && !visitor.time_limit_exceeded(); ++i)
        {
            if (!changed[i])
                continue;

            bool flag = true;
            
            while (flag && !visitor.time_limit_exceeded())
            {
                std::vector<std::size_t> new_set;
                energy_decrease += twocut_lifted::kernighanLin(original_graph, lifted_graph, edgeCosts, partitions[i], new_set, twocut_buffers, twocut_settings);
                
                flag = !new_set.empty();

                if (!new_set.empty())
                    partitions.emplace_back(std::move(new_set));
            }
        }

        if (energy_decrease == .0)
            break;
        
        std::stack<std::size_t> S;
        
        std::fill(visited.begin(), visited.end(), 0);

        // do connected component labeling on the original graph
        numberOfComponents = 0;
        for (std::size_t i = 0; i < original_graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto label = twocut_buffers.vertex_labels[i];

                twocut_buffers.referenced_by[i] = numberOfComponents;

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                        if (twocut_buffers.vertex_labels[it->vertex()] == label && !visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                            twocut_buffers.referenced_by[it->vertex()] = numberOfComponents;
                        }
                }

                ++numberOfComponents;
            }

        // compute new true energy
        double new_energy_value = .0;
        for (std::size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        {
            auto v0 = lifted_graph.vertexOfEdge(i, 0);
            auto v1 = lifted_graph.vertexOfEdge(i, 1);

            if (twocut_buffers.referenced_by[v0] != twocut_buffers.referenced_by[v1])
                new_energy_value += edgeCosts[i];
        }
        
        // if the new true energy is higher, than the current one, revert the changes and terminate
        if (new_energy_value >= current_energy_value - settings.epsilon)
        {
            for (std::size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
            {
                auto v0 = lifted_graph.vertexOfEdge(i, 0);
                auto v1 = lifted_graph.vertexOfEdge(i, 1);

                outputLabels[i] = last_good_vertex_labels[v0] == last_good_vertex_labels[v1] ? 0 : 1;
            }

            break;
        }
        
        // otherwise, form new partitions
        partitions.clear();
        partitions.resize(numberOfComponents);

        for (std::size_t i = 0; i < original_graph.numberOfVertices(); ++i)
        {
            twocut_buffers.vertex_labels[i] = twocut_buffers.referenced_by[i];
            
            partitions[twocut_buffers.vertex_labels[i]].push_back(i);
        }

        twocut_buffers.max_not_used_label = numberOfComponents;

        bool didnt_change = true;
        for (std::size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        {
            auto v0 = lifted_graph.vertexOfEdge(i, 0);
            auto v1 = lifted_graph.vertexOfEdge(i, 1);

            outputLabels[i] = twocut_buffers.vertex_labels[v0] == twocut_buffers.vertex_labels[v1] ? 0 : 1;

            if (static_cast<bool>(outputLabels[i]) != (last_good_vertex_labels[v0] != last_good_vertex_labels[v1]))
                didnt_change = false;
        }

        if (didnt_change)
            break;

        // check if the shape of some partitions didn't change
        changed.clear();
        changed.resize(numberOfComponents);

        std::fill(visited.begin(), visited.end(), 0);

        for (std::size_t i = 0; i < original_graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto label_new = twocut_buffers.vertex_labels[i];
                auto label_old = last_good_vertex_labels[i];

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto w = original_graph.verticesFromVertexBegin(v); w != original_graph.verticesFromVertexEnd(v); ++w)
                    {
                        if (last_good_vertex_labels[*w] == label_old && twocut_buffers.vertex_labels[*w] != label_new)
                            changed[label_new] = 1;

                        if (visited[*w])
                            continue;

                        if (twocut_buffers.vertex_labels[*w] == label_new)
                        {
                            S.push(*w);
                            visited[*w] = 1;

                            if (last_good_vertex_labels[*w] != label_old)
                                changed[label_new] = 1;
                        }
                    }
                }
            }

        last_good_vertex_labels = twocut_buffers.vertex_labels;

        if (!visitor(last_good_vertex_labels) || visitor.time_limit_exceeded())
            break;

        if (settings.verbose)
            std::cout << std::setw(4) << k+1 << std::setw(16) << current_energy_value - new_energy_value << std::setw(15) << energy_decrease << std::setw(15) << ee << std::setw(15) << (energy_decrease - ee) << std::setw(14) << partitions.size() << std::endl;

        current_energy_value = new_energy_value;
    }
}

}
}
}

#endif
