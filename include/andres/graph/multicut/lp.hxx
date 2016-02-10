#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LP_HXX
#define ANDRES_GRAPH_MULTICUT_LP_HXX

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <deque>
#include <algorithm>

#include "andres/graph/complete-graph.hxx"
#include "andres/graph/shortest-paths.hxx"


namespace andres {
namespace graph {
namespace multicut {

template<typename RELAX, typename GRAPH, typename ECA>
inline
std::vector<double> lp(const GRAPH& graph, const ECA& edgeCosts, std::size_t numberOfIterations = std::numeric_limits<std::size_t>::max())
{
    struct Visitor
    {
        bool operator()() const
        {
            return true;
        }
    } visitor;

    return lp<RELAX>(graph, edgeCosts, visitor, numberOfIterations);
}

template<typename RELAX, typename GRAPH, typename ECA, typename VIS>
inline
std::vector<double> lp(const GRAPH& graph, const ECA& edgeCosts, VIS& visitor, std::size_t numberOfIterations = std::numeric_limits<std::size_t>::max())
{
    const double tolerance = std::numeric_limits<float>::epsilon();

    RELAX lp;

    std::deque<std::size_t> path;
    std::vector<double> vars(graph.numberOfEdges());
    std::vector<double> variables(graph.numberOfEdges());
    std::vector<double> coefficients(graph.numberOfEdges());
    std::vector<double> distances(graph.numberOfVertices()); 
    std::vector<std::size_t> parents(graph.numberOfVertices());

    auto addCycleInequalities = [&] ()
    {
        for (std::size_t i = 0; i < vars.size(); ++i)
            vars[i] = std::max(.0, lp.variableValue(i));
            // although in Gurobi we constrain the variables to be in the [0,1] range,
            // sometimes Gurobi finds sligthly negative solutions of the order of 1e-13.
            // The latter totally screws up Dijkstra's shortest path algorithm

        std::size_t counter = 0;

        for (ptrdiff_t edge = 0; edge < graph.numberOfEdges(); ++edge) 
        {
            auto v0 = graph.vertexOfEdge(edge, 0);
            auto v1 = graph.vertexOfEdge(edge, 1);

            // search for shortest path
            double distance;
            spsp(graph, DefaultSubgraphMask<>(), v0, v1, vars.begin(), path, distance, distances.begin(), parents.begin());

            if (vars[edge] > distance + tolerance)
            {
                // add inequality
                auto sz = path.size();

                for (std::size_t j = 0; j < sz - 1; ++j)
                {
                    variables[j] = static_cast<double>(graph.findEdge(path[j], path[j + 1]).second);
                    coefficients[j] = 1.0;
                }

                variables[sz-1] = static_cast<double>(edge);
                coefficients[sz-1] = -1.0;

                lp.addConstraint(variables.begin(), variables.begin() + sz, coefficients.begin(), .0, std::numeric_limits<double>::infinity());

                ++counter;
            }
        }

        return counter;
    };

    lp.initModel(graph.numberOfEdges(), edgeCosts.data());

    for (std::size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        if (!visitor())
            break;

        lp.optimize();

        if (addCycleInequalities() == 0)
            break;
    }

    std::vector<double> edge_values(graph.numberOfEdges());

    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
        edge_values[i] = lp.variableValue(i);

    return edge_values;
}

}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LP_HXX
