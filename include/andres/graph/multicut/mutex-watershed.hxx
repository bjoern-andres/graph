#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_MUTEX_WATERSHED_HXX
#define ANDRES_GRAPH_MULTICUT_MUTEX_WATERSHED_HXX

#include "andres/partition.hxx"

namespace andres {
namespace graph {
namespace multicut {

// Construct a heuristic solution for the Minimum Cost Multicut Problem
// by the mutex watershed algorithm [1].
//
// [1] "The Mutex Watershed: Efficient, Parameter-Free Image Partitioning"
//     S. Wolf, C. Pape, A. Bailoni, N. Rahaman, A. Kreshuk, U. Kothe, F. A. Hamprecht
//     The European Conference on Computer Vision (ECCV), 2018
//    
// Modified implementation by Jan-Hendrik Lange (c) 2018
//
template<typename GRAPH, typename ECA, typename ELA>
void mutexWatershed(GRAPH const& graph, ECA const& edge_costs, ELA& output_labels)
{
    // sort edges for decreasing absolute costs
    std::vector<size_t> edges(graph.numberOfEdges());
    std::iota(edges.begin(), edges.end(), 0);
    std::sort(edges.begin(), edges.end(), [&](const size_t a, const size_t b){
        return std::abs(edge_costs[a]) > std::abs(edge_costs[b]);
    });

    std::vector<std::set<size_t>> mutexes(graph.numberOfVertices());
    andres::Partition<size_t> partition(graph.numberOfVertices());

    /// MAIN LOOP
    for (auto const e : edges)
    {
        auto u = graph.vertexOfEdge(e, 0);
        auto v = graph.vertexOfEdge(e, 1);

        auto ru = partition.find(u);
        auto rv = partition.find(v);

        // if nodes are in the same component
        if (ru == rv)
            continue;

        // if mutex constraint is active
        if (mutexes[ru].find(rv) != mutexes[ru].end())
            continue;

        // otherwise insert mutex constraint for negative edge
        if (edge_costs[e] < 0)
        {
            mutexes[ru].insert(rv);
            mutexes[rv].insert(ru);
        }
        // otherwise merge components for positive edge
        else
        {
            partition.merge(u, v);
            auto ruv = partition.find(ru);
            // update mutex representatives
            if (ruv != ru)
            {
                mutexes[ruv].insert(mutexes[ru].begin(), mutexes[ru].end());
                for (auto & m : mutexes[ru])
                {
                    mutexes[m].erase(ru);
                    mutexes[m].insert(ruv);
                }  
            }
            if (ruv != rv)
            {
                mutexes[ruv].insert(mutexes[rv].begin(), mutexes[rv].end());
                for (auto & m : mutexes[rv])
                {
                    mutexes[m].erase(rv);
                    mutexes[m].insert(ruv);
                }  
            }
        }
    }

    // determine output labels
    for (size_t e = 0; e < graph.numberOfEdges(); e++)
        output_labels[e] = partition.find(graph.vertexOfEdge(e, 0)) != partition.find(graph.vertexOfEdge(e, 1));
}

} // namespace multicut
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_MUTEX_WATERSHED_HXX