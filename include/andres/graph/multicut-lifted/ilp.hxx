#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_ILP_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_ILP_HXX

#include <vector>
#include <deque>
#include <stack>
#include <algorithm>

#include <andres/graph/paths.hxx>
#include <andres/graph/components.hxx>
#include <andres/graph/shortest-paths.hxx>
#include <levinkov/timer.hxx>


namespace andres {
namespace graph {
namespace multicut_lifted
{

template<typename ILP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename ELA>
inline
void ilp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct Visitor
    {
        bool operator()(ELA const& edge_labels) const
        {
            return true;
        }
    } visitor;

    ilp<ILP>(original_graph, lifted_graph, edgeCosts, inputLabels, outputLabels, visitor, numberOfIterations);
}

template<typename ILP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename ELA, typename VIS>
inline
void ilp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, VIS& visitor, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct SubgraphWithCut
    {
        SubgraphWithCut(ILP const& ilp, std::vector<size_t> const& edge_in_lifted_graph)
            : ilp_(ilp), edge_in_lifted_graph_(edge_in_lifted_graph)
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
    levinkov::Timer t;

    std::vector<double> coefficients(lifted_graph.numberOfEdges());
    ComponentsBySearch<ORIGGRAPH> components;
    std::vector<ptrdiff_t> buffer;
    std::vector<size_t> edge_in_lifted_graph(original_graph.numberOfEdges());
    std::deque<size_t> path;
    std::vector<size_t> variables(lifted_graph.numberOfEdges());
    std::vector<char> visited(original_graph.numberOfVertices());

    auto addCycleInequalities = [&] ()
    {
        t.stop();

        levinkov::Timer t_separation;
        t_separation.start();

        components.build(original_graph, SubgraphWithCut(ilp, edge_in_lifted_graph));

        // search for violated non-chordal cycles and add corresp. inequalities
        size_t nCycle = 0;
        size_t nPath = 0;
        size_t nCut = 0;

        for (ptrdiff_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
        {
            auto lv0 = lifted_graph.vertexOfEdge(edge, 0);
            auto lv1 = lifted_graph.vertexOfEdge(edge, 1);

            if (ilp.label(edge) == 1 && components.areConnected(lv0, lv1))
            {
                // if cycle/path inequality is violated

                // search for shortest path
                spsp(original_graph, SubgraphWithCut(ilp, edge_in_lifted_graph), lv0, lv1, path, buffer);

                // check for chords only if the edge is from the original graph
                if (original_graph.findEdge(lv0, lv1).first && findChord(original_graph, path.begin(), path.end(), true).first)
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
            else if (ilp.label(edge) == 0 && !components.areConnected(lv0, lv1))
            {
                // if cut inequality is violated

                // do a simple DFS to find all the edges, that leave the partition, which is essentially a cut
                auto label = components.labels_[lv0];

                std::fill(visited.begin(), visited.end(), 0);
                std::stack<size_t> S;
                S.push(lv0);
                visited[lv0] = 1;

                size_t sz = 0;

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                        if (components.labels_[it->vertex()] == label)
                        {
                            if (!visited[it->vertex()])
                            {
                                S.push(it->vertex());
                                visited[it->vertex()] = 1;
                            }
                        }
                        else
                        {
                            variables[sz] = edge_in_lifted_graph[it->edge()];
                            coefficients[sz++] = -1;
                        }
                }

                // variables[sz] = static_cast<double>(edge);
                // coefficients[sz] = 1;

                // #pragma omp critical
                // ilp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                // #pragma omp atomic
                // ++counter;


                label = components.labels_[lv1];

                S.push(lv1);
                visited[lv1] = 1;

                // sz = 0;

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                        if (components.labels_[it->vertex()] == label)
                        {
                            if (!visited[it->vertex()])
                            {
                                S.push(it->vertex());
                                visited[it->vertex()] = 1;
                            }
                        }
                        else
                        {
                            variables[sz] = edge_in_lifted_graph[it->edge()];
                            coefficients[sz++] = -1;
                        }
                }

                coefficients[sz] = 1;
                variables[sz] = edge;

                ilp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                ++nCut;
            }
        }

        t_separation.stop();

        double objValue = .0;
        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
            objValue += ilp.label(i)*edgeCosts[i];

        std::cerr << t.get_elapsed_seconds() << " " << objValue << " " << nCycle << " " << nPath << " " << nCut << " " << t_separation.get_elapsed_seconds() << std::endl;

        t.start();

        return nCycle + nPath + nCut;
    };

    auto repairSolution = [&] ()
    {
        for(size_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
        {
            auto v0 = lifted_graph.vertexOfEdge(edge, 0);
            auto v1 = lifted_graph.vertexOfEdge(edge, 1);

            outputLabels[edge] = components.areConnected(v0, v1) ? 0 : 1;
        }

        ilp.setStart(outputLabels.begin());
    };

    t.start();

    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto v0 = original_graph.vertexOfEdge(i, 0);
        auto v1 = original_graph.vertexOfEdge(i, 1);

        edge_in_lifted_graph[i] = lifted_graph.findEdge(v0, v1).second;
    }

    ilp.initModel(lifted_graph.numberOfEdges(), edgeCosts.data());
    ilp.setStart(inputLabels.begin());

    for(size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        if (i != 0)
        {
            repairSolution();

            auto energy_value = inner_product(edgeCosts.begin(), edgeCosts.end(), inputLabels.begin(), .0);
            printf("Energy: %f\n", energy_value);
            energy_value = .0;
            for(ptrdiff_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
                if (ilp.label(edge) > 0.5)
                    energy_value += edgeCosts[edge];

            printf("ILP Energy: %f\n", energy_value);

            if (!visitor(outputLabels))
                break;
        }

        ilp.optimize();

        if (addCycleInequalities() == 0)
            break;
    }

    repairSolution();
}

}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_ILP_HXX
