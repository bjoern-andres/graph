#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_ILP_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_ILP_HXX

#include <vector>
#include <deque>
#include <stack>
#include <algorithm>
#include <iomanip>

#include <andres/graph/paths.hxx>
#include <andres/graph/components.hxx>
#include <andres/graph/shortest-paths.hxx>


namespace andres {
namespace graph {
namespace multicut_lifted
{

template<typename ILP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename ELA>
inline
void ilp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct EmptyVisitor
    {
    } visitor;

    ilp<ILP>(original_graph, lifted_graph, edgeCosts, inputLabels, outputLabels, visitor, numberOfIterations);
}

template<typename ILP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename ELA, typename VIS>
inline
void ilp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, VIS& visitor, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct SubgraphWithCut
    {
        SubgraphWithCut(ILP const& ilp, std::vector<size_t> const& edge_in_lifted_graph) :
            ilp_(ilp), edge_in_lifted_graph_(edge_in_lifted_graph)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return ilp_.label(edge_in_lifted_graph_[e]) < .5;
        }

        ILP const& ilp_;
        std::vector<size_t> const& edge_in_lifted_graph_;
    };

    ILP ilp;

    std::vector<double> coefficients(lifted_graph.numberOfEdges());
    ComponentsBySearch<ORIGGRAPH> components;
    std::vector<ptrdiff_t> buffer;
    std::vector<size_t> edge_in_lifted_graph(original_graph.numberOfEdges());
    std::deque<size_t> path;
    std::vector<size_t> variables(lifted_graph.numberOfEdges());
    std::vector<char> visited(original_graph.numberOfVertices());

    auto addCycleInequalities = [&] ()
    {
        components.build(original_graph, SubgraphWithCut(ilp, edge_in_lifted_graph));

        // search for violated non-chordal cycles and add corresp. inequalities
        size_t nCycle = 0;
        size_t nPath = 0;
        size_t nCut = 0;

        for (ptrdiff_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
        {
            auto const lv0 = lifted_graph.vertexOfEdge(edge, 0);
            auto const lv1 = lifted_graph.vertexOfEdge(edge, 1);

            if (ilp.label(edge) > .5 && components.areConnected(lv0, lv1))
            {
                // if cycle/path inequality is violated

                // search for shortest path that contains only non-lifted edges
                spsp(original_graph, SubgraphWithCut(ilp, edge_in_lifted_graph), lv0, lv1, path, buffer);

                bool chordless = true;
                for (auto it1 = path.begin(); it1 != path.end() - 2 && chordless; ++it1)
                    for (auto it2 = it1 + 2; it2 != path.end(); ++it2)
                    {
                        if (it1 == path.begin() && it2 == path.end() - 1)
                            continue;

                        auto const e = lifted_graph.findEdge(*it1, *it2);
                        if (e.first && ilp.label(e.second) > .5)
                        {
                            chordless = false;
                            break;
                        }
                    }

                if (!chordless)
                    continue;

                // add inequality
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    variables[j] = lifted_graph.findEdge(path[j], path[j + 1]).second;
                    coefficients[j] = 1.0;
                }

                variables[path.size() - 1] = edge;
                coefficients[path.size() - 1] = -1.0;

                ilp.addConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), 0, std::numeric_limits<double>::infinity());

                if (original_graph.findEdge(lv0, lv1).first)
                    ++nCycle;
                else
                    ++nPath;
            }
            else if (ilp.label(edge) < .5 && !components.areConnected(lv0, lv1))
            {
                // if cut inequality is violated

                // do a simple DFS to find all the edges, that leave the partition, which is essentially a cut
                auto label = components.labels_[lv0];

                std::fill(visited.begin(), visited.end(), 0);

                std::stack<size_t> S;
                S.push(lv0);
                visited[lv0] = 1;

                ptrdiff_t sz = 0;
                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                        if (components.labels_[it->vertex()] != label)
                        {
                            coefficients[sz] = -1;
                            variables[sz] = edge_in_lifted_graph[it->edge()];
                            
                            ++sz;
                        }
                        else if (!visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                        }
                }

                coefficients[sz] = 1;
                variables[sz] = edge;

                ilp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                ++nCut;


                label = components.labels_[lv1];

                S.push(lv1);
                visited[lv1] = 1;

                sz = 0;
                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                        if (components.labels_[it->vertex()] != label)
                        {
                            coefficients[sz] = -1;
                            variables[sz] = edge_in_lifted_graph[it->edge()];
                            
                            ++sz;
                        }
                        else if (!visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                        }
                }

                coefficients[sz] = 1;
                variables[sz] = edge;

                ilp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                ++nCut;
            }
        }

        double objValue = .0;
        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
            objValue += ilp.label(i)*edgeCosts[i];

        std::cerr << std::fixed << std::setprecision(10) << objValue << " " << nCycle << " " << nPath << " " << nCut << std::endl;

        return nCycle + nPath + nCut;
    };

    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto const v0 = original_graph.vertexOfEdge(i, 0);
        auto const v1 = original_graph.vertexOfEdge(i, 1);

        edge_in_lifted_graph[i] = lifted_graph.findEdge(v0, v1).second;
    }

    ilp.initModel(lifted_graph.numberOfEdges(), edgeCosts.data());
    ilp.setStart(inputLabels.begin());

    for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        ilp.optimize();

        if (addCycleInequalities() == 0)
            break;
    }

    components.build(original_graph, SubgraphWithCut(ilp, edge_in_lifted_graph));

    for (size_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
    {
        auto const v0 = lifted_graph.vertexOfEdge(edge, 0);
        auto const v1 = lifted_graph.vertexOfEdge(edge, 1);

        outputLabels[edge] = components.areConnected(v0, v1) ? 0 : 1;
    }
}

}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_ILP_HXX
