#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_PREPROCESSING_HXX
#define ANDRES_GRAPH_MULTICUT_PREPROCESSING_HXX

#include "andres/graph/components.hxx"
#include "andres/graph/bridges.hxx"

namespace andres {
namespace graph {
namespace multicut {


// Edge contraction operation by masking
// First vertex of given node pair stays in the graph, the other one is removed (masked)
// Flags affected vertices in optional flag pointer
void contract(andres::graph::Graph<> & graph, std::vector<double> & edge_costs, 
    std::vector<char> & emask, std::vector<char> & vmask, std::vector<size_t> pair,
    std::vector<std::pair<std::pair<size_t,size_t>, char>> & constr, std::vector<char> * flag = NULL)
{
    double const tolerance = std::numeric_limits<double>::epsilon();
    auto u = pair[0];
    auto v = pair[1];
    if (!graph.findEdge(u,v).first)
    {
        std::cout << "ERROR (edge contraction): Given pair of vertices is not an edge." << std::endl;
        return;
    }
    auto edge = graph.findEdge(u, v).second;
    if (emask[edge])
    {
        std::cout << "ERROR (edge contraction): Given edge is masked." << std::endl;
        return;
    }
    if (edge_costs[edge] < 0)
    {
        std::cout << "ERROR (edge contraction): Cannot contract, given edge is not positive (" << edge_costs[edge] << ")." << std::endl;
        return;
    }
    emask[edge] = 1;
    if (flag)
        flag->at(u) = 1;
    // update edge costs and remove edges
    for (auto it = graph.adjacenciesFromVertexBegin(v); it != graph.adjacenciesFromVertexEnd(v); it++)
    {
        auto e = it->edge();
        if (emask[e])
            continue;
        auto w = it->vertex();
        if (flag)
            flag->at(w) = 1;
        size_t f;
        if (graph.findEdge(u,w).first)
        {
            f = graph.findEdge(u,w).second;
            // make sure edge is not masked
            if (emask[f])
            {
                std::cout << "WARNING (edge contraction): Edge to be merged was masked." << std::endl;
                emask[f] = 0;
            }
        }
        else
        {
            f = graph.insertEdge(u,w);
            edge_costs.push_back(0.0);
            emask.push_back(0);
        }
        edge_costs[f] += edge_costs[e];
        emask[e] = 1;
    }
    vmask[v] = 1;
}


// Create reduced graph that omits removed edges/vertices
void createReducedGraph(andres::graph::Graph<> & graph, std::vector<double> & edge_costs,
 std::vector<char> & emask, std::vector<char> & vmask, std::vector<size_t> & vorig)
{
    // reserve memory for reduced graph
    andres::graph::Graph<> graph_reduced;
    graph_reduced.reserveVertices(graph.numberOfVertices());
    graph_reduced.reserveEdges(graph.numberOfEdges());
    std::vector<double> edge_costs_reduced;
    edge_costs_reduced.reserve(graph.numberOfEdges());

    std::vector<size_t> vertices_reduced(graph.numberOfVertices());
    std::vector<size_t> vorig_reduced;
    vorig_reduced.reserve(graph.numberOfVertices());

    double const tolerance = std::numeric_limits<double>::epsilon();

    for (size_t v = 0; v < graph.numberOfVertices(); v++)
    {
        if (!vmask[v])
        {
            vertices_reduced[v] = graph_reduced.insertVertex();
            vorig_reduced.push_back(vorig[v]);
        }
    }
    for (size_t e = 0; e < graph.numberOfEdges(); e++)
    {
        if (emask[e] || std::abs(edge_costs[e]) < tolerance)
            continue;
        if (vmask[graph.vertexOfEdge(e, 0)] || vmask[graph.vertexOfEdge(e, 1)])
            continue;
        auto v0 = vertices_reduced[graph.vertexOfEdge(e, 0)];
        auto v1 = vertices_reduced[graph.vertexOfEdge(e, 1)];
        auto e_new = graph_reduced.insertEdge(v0,v1);
        edge_costs_reduced.push_back(edge_costs[e]);

        // perform checks
        assert(v0 < graph_reduced.numberOfVertices());
        assert(v1 < graph_reduced.numberOfVertices());
        assert(e_new == graph_reduced.size()-1);
    }

    // delete isolated vertices
    for (size_t v = 0; v < graph_reduced.numberOfVertices(); v++)
        if (graph_reduced.numberOfEdgesFromVertex(v) == 0)
        {
            vorig_reduced[v] = vorig_reduced[vorig_reduced.size()-1];
            vorig_reduced.pop_back(); 
            graph_reduced.eraseVertex(v);
            v--;
        }

    // set return values
    graph = graph_reduced;
    edge_costs = edge_costs_reduced;
    vorig = vorig_reduced;
}

inline void printRelativeGraphSize(size_t nVertices, size_t nEdges, size_t nOrigVertices, size_t nOrigEdges)
{
    std::cout << "|V| = " << static_cast<double>(nVertices) / static_cast<double>(nOrigVertices) * 100 << "\%"
            << ", |E| = " << static_cast<double>(nEdges) / static_cast<double>(nOrigEdges) * 100 << "\%" << std::endl;
}

// Remove all repulsive edges between components G^+ = (V, E^+)
double extractIsolatedComponents(andres::graph::Graph<> & graph, std::vector<double> & edge_costs,
    std::vector<std::pair<std::pair<size_t,size_t>, char>> & constr, std::vector<size_t> & vorig)
{
    struct AttractionSubgraph
    {
        AttractionSubgraph(std::vector<double> const& edge_costs) :
            edge_costs_(edge_costs)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return (edge_costs_[e] > 0);
        }

        std::vector<double> const& edge_costs_;
    };

    // extract connected components of G^+
    ComponentsBySearch<andres::graph::Graph<>> components;
    components.build(graph, AttractionSubgraph(edge_costs));

    // create reduced graph
    andres::graph::Graph<> graph_reduced;
    graph_reduced.insertVertices(graph.numberOfVertices());
    graph_reduced.reserveEdges(graph.numberOfEdges());
    std::vector<double> edge_costs_reduced;
    edge_costs_reduced.reserve(graph.numberOfEdges());
    
    double cost_offset = 0;
    double const tolerance = std::numeric_limits<double>::epsilon();

    // insert edges
    for (size_t e = 0; e < graph.numberOfEdges(); e++)
    {
        auto u = graph.vertexOfEdge(e, 0);
        auto v = graph.vertexOfEdge(e, 1);

        if (!components.areConnected(u,v))
        {
            // set variable to 1 and remove edge
            cost_offset += edge_costs[e];
            constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 1);
            continue;
        }
        if (std::abs(edge_costs[e]) < tolerance)
            continue;
        auto e_new = graph_reduced.insertEdge(u,v);
        edge_costs_reduced.push_back(edge_costs[e]);
    }

    // delete isolated vertices
    for (size_t v = 0; v < graph_reduced.numberOfVertices(); v++)
        if (graph_reduced.numberOfEdgesFromVertex(v) == 0)
        {
            vorig[v] = vorig[vorig.size()-1];
            vorig.pop_back();
            graph_reduced.eraseVertex(v);
            v--;
        }

    // set return values
    graph = graph_reduced;
    edge_costs = edge_costs_reduced;

    return cost_offset;
}


// Remove all bridges from the graph.
double removeBridges(andres::graph::Graph<> & graph, std::vector<double> & edge_costs,
    std::vector<std::pair<std::pair<size_t,size_t>, char>> & constr, std::vector<size_t> & vorig)
{
    std::vector<char> is_bridge(graph.numberOfEdges());
    double cost_offset = 0;
    findBridges(graph, is_bridge);
    
    // create reduced graph
    andres::graph::Graph<> graph_reduced;
    graph_reduced.insertVertices(graph.numberOfVertices());
    graph_reduced.reserveEdges(graph.numberOfEdges());
    std::vector<double> edge_costs_reduced;
    edge_costs_reduced.reserve(graph.numberOfEdges());

    double const tolerance = std::numeric_limits<double>::epsilon();

    // insert edges
    for (size_t e = 0; e < graph.numberOfEdges(); e++)
    {
        auto u = graph.vertexOfEdge(e, 0);
        auto v = graph.vertexOfEdge(e, 1);

        if (is_bridge[e])
        {
            if (edge_costs[e] < 0)
            {
                cost_offset += edge_costs[e];
                constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 1);
            }
            else
                constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 0);
            continue;
        }
        if (std::abs(edge_costs[e]) < tolerance)
            continue;
        auto e_new = graph_reduced.insertEdge(u,v);
        edge_costs_reduced.push_back(edge_costs[e]);
    }

    // delete isolated vertices
    for (size_t v = 0; v < graph_reduced.numberOfVertices(); v++)
        if (graph_reduced.numberOfEdgesFromVertex(v) == 0)
        {
            vorig[v] = vorig[vorig.size()-1];
            vorig.pop_back();
            graph_reduced.eraseVertex(v);
            v--;
        }

    // set return values
    graph = graph_reduced;
    edge_costs = edge_costs_reduced;

    return cost_offset; 
}


// Repeatedly check single node cuts for dominant edges
double checkSingleNodeCuts(andres::graph::Graph<> & graph, std::vector<double> & edge_costs,
    std::vector<std::pair<std::pair<size_t,size_t>, char>> & constr, std::vector<size_t> & vorig)
{
    std::vector<char> emask(graph.numberOfEdges());
    std::vector<char> vmask(graph.numberOfVertices());
    double cost_offset = 0;
    bool repeat = true;

    std::vector<char> flag(graph.numberOfVertices(), 1);

    while (repeat)
    {
        repeat = false;

        // degree and dominance checks
        for (size_t v = 0; v < graph.numberOfVertices(); v++)
        {
            // only check vertices that exist and are flagged
            if (vmask[v] || !flag[v])
                continue;

            flag[v] = 0;

            // check dominance of incident edges
            double sum_neg_costs = 0;
            double sum_pos_costs = 0;
            double max_abs_cost = -std::numeric_limits<double>::infinity();
            size_t max_edge;
            size_t degree = 0;

            for (auto it = graph.adjacenciesFromVertexBegin(v); it < graph.adjacenciesFromVertexEnd(v); it++)
            {
                auto e = it->edge();
                if (emask[e])
                    continue;
                degree++;

                if (edge_costs[e] > 0)
                    sum_pos_costs += edge_costs[e];
                else
                    sum_neg_costs += edge_costs[e];
                if (std::abs(edge_costs[e]) > max_abs_cost)
                {
                    max_abs_cost = std::abs(edge_costs[e]);
                    max_edge = e;
                }
            }

            // ignore isolated vertices
            if (degree == 0)
            {
                vmask[v] = 1;
                continue;
            }

            // get node u such that {v, u} is max_edge
            size_t u = (graph.vertexOfEdge(max_edge, 0) == v) ? graph.vertexOfEdge(max_edge, 1) : graph.vertexOfEdge(max_edge, 0);

            // check if repulsive edge is dominant
            if (edge_costs[max_edge] < 0 && max_abs_cost >= sum_pos_costs)
            {
                // check if degree 2 vertex
                if (degree == 2)
                {
                    // switch encoding and perform contraction if applicable
                    size_t min_edge;
                    size_t w;
                    // get other edge
                    for (auto it = graph.adjacenciesFromVertexBegin(v); it < graph.adjacenciesFromVertexEnd(v); it++)
                    {
                        if (emask[it->edge()])
                            continue;
                        if (it->edge() != max_edge)
                        {
                            min_edge = it->edge();
                            w = it->vertex();
                        }
                    }
                    // check if other edge is attractive
                    if (edge_costs[min_edge] > 0)
                    {
                        // switch encoding and contract
                        cost_offset += edge_costs[min_edge] + edge_costs[max_edge];
                        edge_costs[min_edge] = -edge_costs[min_edge];
                        edge_costs[max_edge] = -edge_costs[max_edge];
                        contract(graph, edge_costs, emask, vmask, std::vector<size_t>{u, v}, constr, &flag);
                        repeat = true;
                        constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 1);
                    }
                    // otherwise cut and remove both edges
                    else
                    {
                        cost_offset += edge_costs[max_edge] + edge_costs[min_edge];
                        emask[max_edge] = 1;
                        emask[min_edge] = 1;
                        repeat = true;
                        flag[u] = 1;
                        flag[w] = 1;
                        constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 1);
                        constr.emplace_back(std::make_pair(vorig[w],vorig[v]), 1);
                    }
                }
                // cut off degree 1 vertex
                else if (degree == 1)
                {
                    cost_offset += edge_costs[max_edge];
                    emask[max_edge] = 1;
                    repeat = true;
                    flag[u] = 1;
                    constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 1);
                }
            }
            // check if attractive edge is dominant
            else if (edge_costs[max_edge] > 0 && max_abs_cost >= sum_pos_costs - sum_neg_costs - max_abs_cost)
            {
                contract(graph, edge_costs, emask, vmask, std::vector<size_t>{u, v}, constr, &flag);
                repeat = true;
                constr.emplace_back(std::make_pair(vorig[u],vorig[v]), 0);
            }
        }   // end for loop
    }   // end while

    createReducedGraph(graph, edge_costs, emask, vmask, vorig);

    return cost_offset;  
}


// Copyright (c) Jan-Hendrik Lange 2018
//
// Preprocessing routine to shrink an instance of the minimum cost multicut problem.
//
// INPUT:   - graph_orig : Graph of the original instance
//          - edgeCosts  : Edge costs of the original instance
//
// OUTPUT:  - graph : Graph of the reduced instance
//          - edge_costs : Edge costs of the reduced instance
//          - total_cost_offset : Cost of the partial solution
//            (total_cost_offset + cost of solution of reduced instance = cost of corresponding solution of original instance)
//          - constr : Vector of constraints that encode the partial solution.
//              Each entry is a node pair together with a code:
//              { {u,v}, 1} means the original nodes u and v are in distinct clusters
//              { {u,v}, 0} means the original nodes u and v are in the same cluster
//              Note that the nodes u and v need not be adjacent in the original
//              graph. A partial solution is any multicut solution on the original graph that agrees with the constraints.
//          - vorig : Vector that stores for each vertex of the reduced graph the
//                    index of an associated vertex in the original instance
//
template<typename GRAPH, typename ECA>
std::tuple<andres::graph::Graph<>, std::vector<double>, double, std::vector<std::pair<std::pair<size_t, size_t>, char>>, std::vector<size_t>>
preprocessing(GRAPH const& graph_orig, ECA const& edgeCosts, size_t timeLimitSeconds = 86400)
{
    size_t nOrigVertices = graph_orig.numberOfVertices();
    size_t nOrigEdges = graph_orig.numberOfEdges();

    // copy graph
    andres::graph::Graph<> graph;
    graph.insertVertices(graph_orig.numberOfVertices());
    graph.reserveEdges(graph_orig.numberOfEdges());
    for (size_t e = 0; e < graph_orig.numberOfEdges(); e++)
    {
        auto v0 = graph_orig.vertexOfEdge(e, 0);
        auto v1 = graph_orig.vertexOfEdge(e, 1);
        graph.insertEdge(v0,v1);
    }
    // copy edge costs
    std::vector<double> edge_costs(edgeCosts);

    // partial solution constraints
    std::vector<std::pair<std::pair<size_t, size_t>, char>> constr;

    // store for each vertex a corresponding index in the original graph
    std::vector<size_t> vorig(graph.numberOfVertices());
    for (size_t v = 0; v < graph.numberOfVertices(); v++)
        vorig[v] = v;

    double total_cost_offset = 0;
    size_t nEdges = graph.numberOfEdges();
    size_t nVertices = graph.numberOfVertices();
    size_t nEdgesPrev = std::numeric_limits<size_t>::max();

    printRelativeGraphSize(nVertices, nEdges, nOrigVertices, nOrigEdges);

    while (nEdges < nEdgesPrev)
    {
        nEdgesPrev = nEdges;
        total_cost_offset += extractIsolatedComponents(graph, edge_costs, constr, vorig);

        nEdges = graph.numberOfEdges();
        nVertices = graph.numberOfVertices();
        printRelativeGraphSize(nVertices, nEdges, nOrigVertices, nOrigEdges);

        total_cost_offset += checkSingleNodeCuts(graph, edge_costs, constr, vorig);

        nEdges = graph.numberOfEdges();
        nVertices = graph.numberOfVertices();
        printRelativeGraphSize(nVertices, nEdges, nOrigVertices, nOrigEdges);
    }

    return std::make_tuple(graph, edge_costs, total_cost_offset, constr, vorig);  
}


}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_PREPROCESSING_HXX
