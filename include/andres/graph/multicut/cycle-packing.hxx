#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_CYCLE_PACKING_HXX
#define ANDRES_GRAPH_MULTICUT_CYCLE_PACKING_HXX

#include <deque>
#include <algorithm>

#include "andres/graph/components.hxx"
#include "andres/graph/shortest-paths.hxx"
#include "andres/graph/graph.hxx"

namespace andres {
namespace graph {
namespace multicut {


// Copyright (c) Jan-Hendrik Lange 2018
//
// Heuristic algorithm that computes a dual lower bound
// by iteratively packing conflicted cycles.
//
// OUTPUT:  - dual lower bound
//          - reduced primal cost vector associated with heuristic dual solution
//
template<typename GRAPH, typename ECA>
std::tuple<double, ECA>
iterativeCyclePacking(GRAPH const& graph_orig, ECA const& edgeCosts, bool verbose = false, size_t max_cycle_length = 0)
{
    std::vector<double> edge_costs(edgeCosts);

    std::vector<std::pair<size_t, size_t>> repulsive_edges;

    double lower_bound = 0;
    double const tolerance = std::numeric_limits<double>::epsilon();

    // copy attraction subgraph
    andres::graph::Graph<> graph(graph_orig.numberOfVertices());
    graph.reserveEdges(graph_orig.numberOfEdges());

    for (size_t e = 0; e < graph_orig.numberOfEdges(); e++)
    {
        auto v0 = graph_orig.vertexOfEdge(e, 0);
        auto v1 = graph_orig.vertexOfEdge(e, 1);

        if (edge_costs[e] > tolerance)
            graph.insertEdge(v0,v1);
        if (edge_costs[e] < -tolerance)
        {
            repulsive_edges.emplace_back(v0, v1);
            lower_bound += edge_costs[e];
        }
    }
    
    if (verbose)
    {
        std::cout << "Conflicted cycle packing.." << std::endl;
        double trivial_bound = lower_bound;
        std::cout << "Trivial bound: " << trivial_bound << std::endl;
    }

    ComponentsBySearch<andres::graph::Graph<>> components;

    std::deque<size_t> path;
    std::vector<ptrdiff_t> parents(graph.numberOfVertices());
    std::vector<size_t> labels(graph.numberOfVertices());
    size_t seen_label = 0;

    size_t nCycles = 0;

    // check optional cycle length parameter
    if (max_cycle_length < 3)
        max_cycle_length = graph.numberOfEdges();

    // iterative length cycle search
    for (size_t cycle_length = 3; cycle_length <= max_cycle_length; cycle_length++)
    {
        if (verbose)
            std::cout << "Round " << cycle_length-3 << ", L = " << lower_bound << std::endl;

        // update component labeling
        components.build(graph);

        size_t progress = 0;

        // shuffling can give great speed-up if edges with similar indices are spatially close in the graph
        std::random_shuffle(repulsive_edges.begin(), repulsive_edges.end());

        std::vector<std::pair<size_t, size_t>> rem_repuls_edges;
        rem_repuls_edges.reserve(repulsive_edges.size());

        // iterate over repulsive edges
        for (auto p : repulsive_edges)
        {
            // periodically update component labeling for speed-up
            progress++;
            if (progress > 0.1 * repulsive_edges.size())
            {
                components.build(graph);
                progress = 0;
            }

            auto v0 = p.first;
            auto v1 = p.second;
            size_t f = graph_orig.findEdge(v0, v1).second;
            
            // check if conflicted cycle exists
            if (!components.areConnected(v0, v1))
                continue;

            // pack short cycles as long as available and positive weight left
            while (edge_costs[f] < -tolerance)
            {
                seen_label++;
                if (seen_label == 0)
                {
                    std::fill(labels.begin(), labels.end(), 0);
                    seen_label++;
                }
                if (!spsp(graph, DefaultSubgraphMask<>(), v0, v1, path, parents, labels, seen_label, cycle_length-1))
                    break;
                
                // find minimum weight edge in cycle
                double min_weight = -edge_costs[f];
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    auto e = graph_orig.findEdge(path[j], path[j + 1]).second;
                    if (edge_costs[e] < min_weight)
                        min_weight = edge_costs[e];
                }
                // subtract minimum weight and delete edges of negligible weight
                edge_costs[f] += min_weight;
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    auto e = graph_orig.findEdge(path[j], path[j + 1]).second;
                    edge_costs[e] -= min_weight;
                    if (edge_costs[e] < tolerance)
                        graph.eraseEdge(graph.findEdge(path[j], path[j+1]).second);
                }
                // update lower bound
                lower_bound += min_weight;
                nCycles++;             
            }

            if (edge_costs[f] < -tolerance)
                rem_repuls_edges.emplace_back(v0, v1);
        }

        repulsive_edges = rem_repuls_edges;
        if (repulsive_edges.empty())
            break;
    }

    if (verbose)
    {
        std::cout << "Improved bound: " << lower_bound << std::endl;
        std::cout << "Packed " << nCycles << " conflicted cycles." << std::endl;
    }

    return std::make_tuple(lower_bound, edge_costs);    
}


}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_CYCLE_PACKING_HXX
