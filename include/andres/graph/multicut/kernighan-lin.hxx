#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_KERNIGHAN_LIN_HXX
#define ANDRES_GRAPH_MULTICUT_KERNIGHAN_LIN_HXX

#include <iomanip>
#include <stack>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "andres/graph/components.hxx"
#include "andres/graph/twocut/kernighan-lin.hxx"



namespace andres {
namespace graph {
namespace multicut {

struct Settings
{
    std::size_t numberOfInnerIterations { std::numeric_limits<std::size_t>::max() };
    std::size_t numberOfOuterIterations { 100 };
    double epsilon { 1e-7 };
    bool verbose { true };
};

template<typename GRAPH, typename ECA, typename ELA>
inline
void kernighanLin(const GRAPH& graph, const ECA& edgeCosts, const ELA& inputLabels, ELA& outputLabels, const Settings settings = Settings())
{
    struct Visitor
    {
        bool operator()(std::vector<std::size_t>& vertex_labels) const
        {
            return true;
        }
    } visitor;

    kernighanLin(graph, edgeCosts, inputLabels, outputLabels, visitor, settings);
}

template<typename GRAPH, typename ECA, typename ELA, typename VIS>
inline
void kernighanLin(const GRAPH& graph, const ECA& edgeCosts, const ELA& inputLabels, ELA& outputLabels, VIS& visitor, const Settings settings = Settings())
{
    struct SubgraphWithCut { // a subgraph with cut mask
        SubgraphWithCut(const ELA& labels)
            : labels_(labels) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return labels_[e] == 0; }

        const ELA& labels_;
    };

    // build decomposition based on the current multicut
    ComponentsBySearch<GRAPH> components;
    components.build(graph, SubgraphWithCut(inputLabels));

    double starting_energy = .0;

    // find out how many connected components there are
    // check if the input multicut labeling is valid
    std::size_t numberOfComponents = 0;
    for(std::size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
    {
        outputLabels[edge] = inputLabels[edge];

        auto v0 = graph.vertexOfEdge(edge, 0);
        auto v1 = graph.vertexOfEdge(edge, 1);

        numberOfComponents = std::max(numberOfComponents, std::max(components.labels_[v0], components.labels_[v1]));

        if (inputLabels[edge])
            starting_energy += edgeCosts[edge];

        if (static_cast<bool>(inputLabels[edge]) != !components.areConnected(v0, v1))
            throw std::runtime_error("the input multicut labeling is invalid.");
    }

    ++numberOfComponents;

    auto twocut_buffers = twocut::makeTwoCutBuffers(graph);

    twocut::TwoCutSettings twocut_settings;
    twocut_settings.numberOfIterations = settings.numberOfInnerIterations;
    twocut_settings.epsilon = settings.epsilon;

    // build partitions
    std::vector<std::vector<std::size_t>> partitions(numberOfComponents);

    for (std::size_t i = 0; i < components.labels_.size(); ++i)
    {
        partitions[components.labels_[i]].push_back(i);
        twocut_buffers.vertex_labels[i] = components.labels_[i];
    }
    
    twocut_buffers.max_not_used_label = partitions.size();

    if (settings.verbose)
    {
        std::cout << "Starting energy: " << starting_energy << std::endl;
        std::cout << std::setw(4) << "Iter" << std::setw(16) << "Total decrease" << std::setw(15) << "Pair updates" << std::setw(15) << "New sets" << std::setw(15) << "Num. of sets\n";
    }

    auto last_good_vertex_labels = twocut_buffers.vertex_labels;

    // auxillary array for BFS/DFS
    std::vector<char> visited(graph.numberOfVertices());

    // 1 if i-th partitioned changed since last iteration, 0 otherwise
    std::vector<char> changed(numberOfComponents, 1);

    // interatively update bipartition in order to minimize the total cost of the multicut
    for (std::size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {
        auto energy_decrease = .0;

        std::vector<std::unordered_set<std::size_t>> edges(numberOfComponents);
        for (std::size_t e = 0; e < graph.numberOfEdges(); ++e)
            if (outputLabels[e])
            {
                auto v0 = twocut_buffers.vertex_labels[graph.vertexOfEdge(e, 0)];
                auto v1 = twocut_buffers.vertex_labels[graph.vertexOfEdge(e, 1)];

                edges[std::min(v0, v1)].insert(std::max(v0, v1));
            }

        for (std::size_t i = 0; i < numberOfComponents; ++i)
            if (!partitions[i].empty())
                for (auto j : edges[i])
                    if (!partitions[j].empty() && (changed[j] || changed[i]))
                    {
                        auto ret = twocut::kernighanLin(graph, edgeCosts, partitions[i], partitions[j], twocut_buffers, twocut_settings);

                        if (ret > settings.epsilon)
                            changed[i] = changed[j] = 1;

                        energy_decrease += ret;

                        if (partitions[i].size() == 0)
                            break;
                    }
        
        auto ee = energy_decrease;

        // remove partitions that became empty after the previous step
        auto new_end = std::partition(partitions.begin(), partitions.end(), [](const std::vector<std::size_t>& s) { return !s.empty(); });
        partitions.resize(new_end - partitions.begin());

        // try to intoduce new partitions
        for (std::size_t i = 0, p_size = partitions.size(); i < p_size; ++i)
        {
            if (!changed[i])
                continue;

            bool flag = true;
            
            while (flag)
            {
                std::vector<std::size_t> new_set;
                energy_decrease += twocut::kernighanLin(graph, edgeCosts, partitions[i], new_set, twocut_buffers, twocut_settings);

                flag = !new_set.empty();

                if (!new_set.empty())
                    partitions.emplace_back(std::move(new_set));
            }
        }

        if (!visitor(twocut_buffers.vertex_labels))
            break;

        if (energy_decrease == .0)
            break;

        std::stack<std::size_t> S;
        
        std::fill(visited.begin(), visited.end(), 0);

        partitions.clear();
        numberOfComponents = 0;

        // do connected component labeling on the original graph and form new pфrtitions
        for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto label = twocut_buffers.vertex_labels[i];

                twocut_buffers.referenced_by[i] = numberOfComponents;

                partitions.emplace_back(std::vector<std::size_t>());
                partitions.back().push_back(i);

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto it = graph.adjacenciesFromVertexBegin(v); it != graph.adjacenciesFromVertexEnd(v); ++it)
                        if (twocut_buffers.vertex_labels[it->vertex()] == label && !visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                            twocut_buffers.referenced_by[it->vertex()] = numberOfComponents;
                            partitions.back().push_back(it->vertex());
                        }
                }

                ++numberOfComponents;
            }

        twocut_buffers.vertex_labels = twocut_buffers.referenced_by;

        twocut_buffers.max_not_used_label = numberOfComponents;

        bool didnt_change = true;
        for (std::size_t i = 0; i < graph.numberOfEdges(); ++i)
        {
            auto v0 = graph.vertexOfEdge(i, 0);
            auto v1 = graph.vertexOfEdge(i, 1);

            outputLabels[i] = twocut_buffers.vertex_labels[v0] == twocut_buffers.vertex_labels[v1] ? 0 : 1;

            if (static_cast<bool>(outputLabels[i]) != (last_good_vertex_labels[v0] != last_good_vertex_labels[v1]))
                didnt_change = false;
        }

        if (didnt_change)
            break;

        // check if the shape of some partitions didn't change
        changed.resize(numberOfComponents);
        std::fill(changed.begin(), changed.end(), 0);

        std::fill(visited.begin(), visited.end(), 0);

        for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
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

                    for (auto w = graph.verticesFromVertexBegin(v); w != graph.verticesFromVertexEnd(v); ++w)
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

        if (settings.verbose)
            std::cout << std::setw(4) << k+1 << std::setw(16) << energy_decrease << std::setw(15) << ee << std::setw(15) << (energy_decrease - ee) << std::setw(14) << partitions.size() << std::endl;
    }
}

template<typename GraphVisitor, typename ECA, typename ELA>
inline
void kernighanLin(const CompleteGraph<GraphVisitor>& graph, const ECA& edgeCosts, const ELA& inputLabels, ELA& outputLabels, const Settings settings = Settings())
{
    struct SubgraphWithCut { // a subgraph with cut mask
        SubgraphWithCut(const ELA& labels)
            : labels_(labels) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return labels_[e] == 0; }

        const ELA& labels_;
    };

    // build decomposition based on the current multicut
    ComponentsBySearch<CompleteGraph<GraphVisitor>> components;
    components.build(graph, SubgraphWithCut(inputLabels));

    // find out how many connected components there are
    // check if the input multicut labeling is valid
    std::size_t numberOfComponents = 0;
    for(std::size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
    {
        auto v0 = graph.vertexOfEdge(edge, 0);
        auto v1 = graph.vertexOfEdge(edge, 1);

        numberOfComponents = std::max(numberOfComponents, std::max(components.labels_[v0], components.labels_[v1]));

        if (static_cast<bool>(inputLabels[edge]) != !components.areConnected(v0, v1))
            throw std::runtime_error("the input multicut labeling is invalid.");
    }

    ++numberOfComponents;

    std::vector<std::vector<std::size_t>> partitions(numberOfComponents);

    auto twocut_buffers = twocut::makeTwoCutBuffers(graph);
    
    twocut::TwoCutSettings twocut_settings;
    twocut_settings.numberOfIterations = settings.numberOfInnerIterations;
    twocut_settings.epsilon = settings.epsilon;

    // build partitions
    for (std::size_t i = 0; i < components.labels_.size(); ++i)
        partitions[components.labels_[i]].push_back(i);

    // interatively update bipartition in order to minimize the total cost of the multicut
    for (std::size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {
        auto energy_decrease = .0;

        // update pairs of partitions
        for (std::size_t i = 0; i < partitions.size() - 1; ++i)
            for (auto j = i + 1; j < partitions.size(); ++j)
                if (!partitions[j].empty())
                    energy_decrease += twocut::kernighanLin(graph, edgeCosts, partitions[i], partitions[j], twocut_buffers, twocut_settings);

        // remove partitions that became empty after the previous step
        auto new_end = std::partition(partitions.begin(), partitions.end(), [](const std::vector<std::size_t>& s) { return !s.empty(); });
        partitions.resize(new_end - partitions.begin());

        // try to intoduce new partitions
        for (std::size_t i = 0, p_size = partitions.size(); i < p_size; ++i)
        {
            std::vector<std::size_t> new_set;
            energy_decrease += twocut::kernighanLin(graph, edgeCosts, partitions[i], new_set, twocut_buffers, twocut_settings);

            if (!new_set.empty())
                partitions.emplace_back(std::move(new_set));
        }

        if (energy_decrease == .0)
            break;

        for (std::size_t i = 0; i < partitions.size(); ++i)
            for (std::size_t a = 0; a < partitions[i].size(); ++a)
            {
                for (std::size_t b = a + 1; b < partitions[i].size(); ++b)
                    outputLabels[graph.findEdge(partitions[i][a], partitions[i][b]).second] = 0;

                for (std::size_t j = i + 1; j < partitions.size(); ++j)
                    for (auto b : partitions[j])
                        outputLabels[graph.findEdge(partitions[i][a], b).second] = 1;                        
            }
    }
}

} // of multicut
} // of graph
} // of andres

#endif
