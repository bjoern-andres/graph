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
// If oddWheelPacking is true, then the algorithm iteratively
// packs conflicted odd wheels before the default algorithm starts.
//
// OUTPUT:  - dual lower bound
//          - reduced primal cost vector associated with heuristic dual solution
//
template<typename GRAPH, typename ECA>
std::tuple<double, ECA>
iterativeCyclePacking(GRAPH const& graph_orig, ECA const& edgeCosts, bool oddWheelPacking = false, bool verbose = false, size_t max_cycle_length = 0)
{
    
    struct AttractionSubgraph
    {
        AttractionSubgraph(ECA const& edge_costs, std::vector<size_t> & orig_index) :
            edge_costs_(edge_costs), orig_index_(orig_index)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return (edge_costs_[orig_index_[e]] > 0);
        }

        ECA const& edge_costs_;
        std::vector<size_t> & orig_index_;
    };

    // Edge index maintenance
    struct EdgeDeletionVisitor {

        EdgeDeletionVisitor(std::vector<size_t> & orig_index, std::vector<size_t> & new_index, ECA & weights) : 
            orig_index_(orig_index), new_index_(new_index), weights_(weights)
        {}

        void insertVertex(const size_t a) const {}
        void insertVertices(const size_t a, const size_t n) const {}
        void eraseVertex(const size_t a) const {}
        void relabelVertex(const size_t a, const size_t b) const {}
        void insertEdge(const size_t a) const {}
        void eraseEdge(const size_t a) const
        {
            auto last_weight = weights_.size() - 1;
            weights_[a] = weights_[last_weight];
            weights_.pop_back();
        }
        void relabelEdge(const size_t former, const size_t future) const
        {
            orig_index_[future] = orig_index_[former];
            new_index_[orig_index_[future]] = future;
        }

        // stores at position i the original index of edge i
        std::vector<size_t> & orig_index_;
        // stores at position i the new index of original edge i
        std::vector<size_t> & new_index_;
        // weights of dynamic graph
        ECA & weights_;
    };

    std::vector<double> edge_costs(edgeCosts);

    std::set<size_t> repulse_edges;

    double lower_bound = 0;
    double const tolerance = std::numeric_limits<double>::epsilon();

    // create visitor that keeps track of edge indices and weights
    std::vector<size_t> orig_index(graph_orig.numberOfEdges());
    std::vector<size_t> new_index(graph_orig.numberOfEdges());
    std::vector<double> weights(graph_orig.numberOfEdges());

    for (size_t e = 0; e < graph_orig.numberOfEdges(); e++)
    {
        orig_index[e] = e;
        new_index[e] = e;

        // set up weight
        if (edge_costs[e] > tolerance)
        {
            weights[e] = edge_costs[e];

        }
        else if (edge_costs[e] < -tolerance)
        {
            weights[e] = -edge_costs[e];
            repulse_edges.insert(e);
            lower_bound += edge_costs[e];
        }
        else
            weights[e] = 0;
    }

    EdgeDeletionVisitor visitor(orig_index, new_index, weights);
    AttractionSubgraph attr_subgraph(edge_costs, orig_index);

    // copy graph
    andres::graph::Graph<EdgeDeletionVisitor> graph(visitor);
    graph.insertVertices(graph_orig.numberOfVertices());
    for (size_t e = 0; e < graph_orig.numberOfEdges(); e++)
    {
        auto v0 = graph_orig.vertexOfEdge(e, 0);
        auto v1 = graph_orig.vertexOfEdge(e, 1);
        graph.insertEdge(v0,v1);
    }

    /// USAGE OF EDGE INDEX MAINTENANCE STARTS HERE ///

    // delete edges of zero weight
    for (size_t e = 0; e < graph.numberOfEdges(); e++)
    {
        if (weights[e] < tolerance)
        {
            graph.eraseEdge(e);
            e--;
        }
    }
    
    if (verbose)
    {
        std::cout << "Conflicted cycle packing.." << std::endl;
        double trivial_bound = lower_bound;
        std::cout << "Trivial bound: " << trivial_bound << std::endl;
    }

    ComponentsBySearch<andres::graph::Graph<EdgeDeletionVisitor>> components;

    std::deque<size_t> path;
    std::vector<ptrdiff_t> buffer;

    // optional odd-wheel packing
    if (oddWheelPacking)
    {
        size_t nAddWheels = 0;
        std::vector<size_t> vertices(graph.numberOfVertices());
        for (size_t v = 0; v < graph.numberOfVertices(); v++)
            vertices[v] = v;

        std::random_shuffle(vertices.begin(), vertices.end());

        // iterate over vertices
        for (size_t i = 0; i < graph.numberOfVertices(); i++)
        {
            auto v = vertices[i];

            std::vector<char> active_vertex(graph.numberOfVertices());
            std::vector<size_t> neighbors;
            size_t n;

            // get attractive neighborhood
            std::vector<size_t> pred_tree(graph.numberOfVertices());
            std::vector<size_t> distances(graph.numberOfVertices());
            std::fill(distances.begin(), distances.end(), graph.numberOfVertices()+1);

            std::queue<size_t> qtree;
            pred_tree[v] = v;
            distances[v] = 0;
            qtree.push(v);

            size_t max_distance = 1;

            while (!qtree.empty())
            {
                auto n = qtree.front();
                qtree.pop();

                for (auto it = graph.adjacenciesFromVertexBegin(n); it != graph.adjacenciesFromVertexEnd(n); it++)
                {
                    auto e = it->edge();
                    auto w = it->vertex();

                    if (edge_costs[orig_index[e]] < 0)
                        continue;

                    if (distances[w] > graph.numberOfVertices())
                    {
                        distances[w] = distances[n] + 1;
                        pred_tree[w] = n;
                        if (distances[w] < max_distance)
                            qtree.push(w);
                        if (distances[w] == max_distance-1 || distances[w] == max_distance)
                        {
                            neighbors.push_back(w);
                            active_vertex[w] = 1;
                        }
                    }
                }
            }

            // minimal number of neighbors needed for odd cycle
            if (neighbors.size() < 3)
                continue;

            // iterate over neighbors
            for (auto nit = neighbors.begin(); nit != neighbors.end(); nit++)
            {
                n = *nit;
                if (!active_vertex[n])
                    continue;
                bool repeat = false;

                // search for odd cycle in neighbor subgraph
                std::vector<size_t> pred(graph.numberOfVertices());
                std::vector<signed int> color(graph.numberOfVertices());
                std::queue<size_t> Q;
                Q.push(n);
                color[n] = 1;
                pred[n] = n;

                while (!Q.empty() && !repeat)
                {
                    auto u = Q.front();
                    Q.pop();
                    for (auto it = graph.adjacenciesFromVertexBegin(u); it != graph.adjacenciesFromVertexEnd(u); it++)
                    {
                        auto e = it->edge();
                        auto w = it->vertex();

                        if (edge_costs[orig_index[e]] > 0 || weights[e] < tolerance || !active_vertex[w])
                            continue;

                        if (!color[w])
                        {
                            color[w] = -1 * color[u];
                            pred[w] = u;
                            Q.push(w);
                        }

                        // found odd cycle
                        if (color[w] == color[u])
                        {
                            std::vector<std::pair<size_t,size_t>> cycle_edges;
                            std::vector<std::pair<size_t,size_t>> spokes;
                            cycle_edges.emplace_back(graph.vertexOfEdge(e,0), graph.vertexOfEdge(e, 1));
                            double min_weight = weights[e];
                            // traverse back first half of cycle
                            auto p = u;
                            while (p != pred[p])
                            {
                                auto cycle_edge = graph.findEdge(pred[p],p).second;
                                cycle_edges.emplace_back(pred[p], p);
                                if (weights[cycle_edge] < min_weight)
                                    min_weight = weights[cycle_edge];

                                // traverse back leg of star tree
                                size_t s = p;
                                while (s != pred_tree[s])
                                {
                                    auto spoke = graph.findEdge(s,pred_tree[s]).second;
                                    spokes.emplace_back(s, pred_tree[s]);
                                    if (weights[spoke] < min_weight)
                                        min_weight = weights[spoke];
                                    s = pred_tree[s];
                                }

                                p = pred[p];
                            }
                            // traverse back second half of cycle
                            p = w;
                            while (p != pred[p])
                            {
                                auto cycle_edge = graph.findEdge(pred[p],p).second;
                                cycle_edges.emplace_back(pred[p], p);
                                if (weights[cycle_edge] < min_weight)
                                    min_weight = weights[cycle_edge];

                                // traverse back leg of star tree
                                size_t s = p;
                                while (s != pred_tree[s])
                                {
                                    auto spoke = graph.findEdge(s,pred_tree[s]).second;
                                    spokes.emplace_back(s, pred_tree[s]);
                                    if (weights[spoke] < min_weight)
                                        min_weight = weights[spoke];
                                    s = pred_tree[s];
                                }

                                p = pred[p];
                            }
                            // traverse back leg of star tree for start vertex
                            size_t s = p;
                            while (s != pred_tree[s])
                            {
                                auto spoke = graph.findEdge(s,pred_tree[s]).second;
                                spokes.emplace_back(s, pred_tree[s]);
                                if (weights[spoke] < min_weight)
                                    min_weight = weights[spoke];
                                s = pred_tree[s];
                            }

                            // pack odd wheel
                            lower_bound += 0.5 * (static_cast<double>(cycle_edges.size()) + 1.0) * min_weight;

                            for (auto vpair : spokes)
                            {
                                auto edge = graph.findEdge(vpair.first, vpair.second).second;
                                weights[edge] -= min_weight;
                                edge_costs[orig_index[edge]] -= min_weight;
                                if (weights[edge] < tolerance)
                                {
                                    // remove vertices from neighbors
                                    auto v0 = graph.vertexOfEdge(edge, 0);
                                    auto v1 = graph.vertexOfEdge(edge, 1);
                                    active_vertex[v0] = 0;
                                    active_vertex[v1] = 0;
                                    // delete edge
                                    graph.eraseEdge(edge);
                                }

                            }
                            for (auto vpair : cycle_edges)
                            {
                                auto edge = graph.findEdge(vpair.first, vpair.second).second;
                                weights[edge] -= min_weight;
                                edge_costs[orig_index[edge]] += min_weight;
                                if (weights[edge] < tolerance)
                                {
                                    repulse_edges.erase(orig_index[edge]);
                                    graph.eraseEdge(edge);
                                }
                            }
                    
                            nAddWheels++;
                            repeat = true;
                            break;
                        }
                    }
                }

                if (repeat)
                    nit--;
                
            } // end neighbors loop

        } // end vertices loop

        if (verbose)
            std::cout << "Packed " << nAddWheels << " conflicted odd wheels." << std::endl;

    } // end odd wheel packing
    

    // copy remaining repulsive edges to array
    std::vector<size_t> repulsive_edges;
    for (auto it = repulse_edges.begin(); it != repulse_edges.end(); it++)
    {
        repulsive_edges.push_back(*it);
    }

    size_t nCycles = 0;

    // check optional cycle length parameter
    if (max_cycle_length < 3)
        max_cycle_length = graph.numberOfEdges();

    // iterative length cycle search
    for (size_t cycle_length = 3; cycle_length <= max_cycle_length; cycle_length++)
    {
        // shuffling can give great speed-up if edges with similar indices are spatially close in the graph
        std::random_shuffle(repulsive_edges.begin(), repulsive_edges.end());

        if (verbose)
            std::cout << "Round " << cycle_length-3 << ", L = " << lower_bound << std::endl;

        // update component labeling
        components.build(graph, attr_subgraph);

        size_t progress = 0;

        // iterate over repulsive edges
        for (size_t i = 0; i < repulsive_edges.size(); i++)
        {
            // periodically update component labeling for speed-up
            progress++;
            if (progress > 0.05 * repulsive_edges.size())
            {
                components.build(graph, attr_subgraph);
                progress = 0;
            }
        
            auto f = new_index[repulsive_edges[i]];
            auto v0 = graph.vertexOfEdge(f, 0);
            auto v1 = graph.vertexOfEdge(f, 1);
            
            // check if conflicted cycle exists and repulsive edge has positive weight
            if (!components.areConnected(v0, v1) || weights[f] < tolerance)
            {
                auto imax = repulsive_edges.size() - 1;
                repulsive_edges[i] = repulsive_edges[imax];
                repulsive_edges.pop_back();
                i--;
                if (graph.findEdge(v0, v1).first)
                    graph.eraseEdge(f);
                continue;
            }

            // balance short cycles as long as available and positive weight left
            while (weights[f] >= tolerance && spsp(graph, attr_subgraph, v0, v1, path, buffer, cycle_length-1))
            {
                nCycles++;
                // find minimum weight edge in cycle
                double min_weight = weights[f];
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    auto e = graph.findEdge(path[j], path[j + 1]).second;
                    if (weights[e] < min_weight)
                        min_weight = weights[e];
                }
                // subtract minimum weight and remove edges of negligible weight
                weights[f] -= min_weight;
                edge_costs[orig_index[f]] += min_weight;
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    auto e = graph.findEdge(path[j], path[j + 1]).second;
                    weights[e] -= min_weight;
                    edge_costs[orig_index[e]] -= min_weight;
                    if (weights[e] < tolerance)
                        graph.eraseEdge(e);
                }
                // update lower bound
                lower_bound += min_weight;
                // remove repulsive edge if weight is negligible
                if (weights[f] < tolerance)
                {
                    auto imax = repulsive_edges.size() - 1;
                    repulsive_edges[i] = repulsive_edges[imax];
                    repulsive_edges.pop_back();
                    i--;
                    graph.eraseEdge(f);
                    break;
                }
                // update index of repulsive edge (it may change after edge deletions)
                f = graph.findEdge(v0, v1).second;                
            }
        }     

        // terminate if no conflicted cycle remains
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
