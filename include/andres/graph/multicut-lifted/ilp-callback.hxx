#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_ILP_CALLBACK_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_ILP_CALLBACK_HXX

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
void ilp_callback(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, size_t timeLimitSeconds = 86400)
{
    struct EmptyVisitor
    {
    } visitor;

    ilp_callback<ILP>(original_graph, lifted_graph, edgeCosts, inputLabels, outputLabels, visitor, timeLimitSeconds);
}

template<typename ILP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename ELA, typename VIS>
inline
void ilp_callback(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, VIS& visitor, size_t timeLimitSeconds = 86400)
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
    
    class Callback: public ILP::Callback
    {
    public:
        Callback(ILP& solver, ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph) :
            ILP::Callback(solver),
            original_graph_(original_graph),
            lifted_graph_(lifted_graph),
            coefficients_(lifted_graph.numberOfEdges()),
            variables_(lifted_graph.numberOfEdges()),
            visited_(original_graph.numberOfVertices()),
            edge_in_lifted_graph_(original_graph.numberOfEdges())
        {
            for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
            {
                auto const v0 = original_graph.vertexOfEdge(i, 0);
                auto const v1 = original_graph.vertexOfEdge(i, 1);

                edge_in_lifted_graph_[i] = lifted_graph.findEdge(v0, v1).second;
            }
        }

        void separateAndAddLazyConstraints() override
        {
            components_.build(original_graph_, SubgraphWithCut(*this, edge_in_lifted_graph_));

            // search for violated non-chordal cycles and add corresp. inequalities
            size_t nCycle = 0;
            size_t nPath = 0;
            size_t nCut = 0;

            for (ptrdiff_t edge = 0; edge < lifted_graph_.numberOfEdges(); ++edge)
            {
                auto const lv0 = lifted_graph_.vertexOfEdge(edge, 0);
                auto const lv1 = lifted_graph_.vertexOfEdge(edge, 1);

                if (this->label(edge) == 1 && components_.areConnected(lv0, lv1))
                {
                    // if cycle/path inequality is violated

                    // search for shortest path that contains only non-lifted edges
                    spsp(original_graph_, SubgraphWithCut(*this, edge_in_lifted_graph_), lv0, lv1, path_, buffer_);

                    bool chordless = true;
                    for (auto it1 = path_.begin(); it1 != path_.end() - 2 && chordless; ++it1)
                        for (auto it2 = it1 + 2; it2 != path_.end(); ++it2)
                        {
                            if (it1 == path_.begin() && it2 == path_.end() - 1)
                                continue;

                            auto const e = lifted_graph_.findEdge(*it1, *it2);
                            if (e.first && this->label(e.second) > .5)
                            {
                                chordless = false;
                                break;
                            }
                        }

                    if (!chordless)
                        continue;

                    // add inequality
                    for (size_t j = 0; j < path_.size() - 1; ++j)
                    {
                        variables_[j] = lifted_graph_.findEdge(path_[j], path_[j + 1]).second;
                        coefficients_[j] = 1.0;
                    }

                    variables_[path_.size() - 1] = edge;
                    coefficients_[path_.size() - 1] = -1.0;

                    this->addLazyConstraint(variables_.begin(), variables_.begin() + path_.size(), coefficients_.begin(), 0, std::numeric_limits<double>::infinity());

                    if (original_graph_.findEdge(lv0, lv1).first)
                        ++nCycle;
                    else
                        ++nPath;
                }
                else if (this->label(edge) == 0 && !components_.areConnected(lv0, lv1))
                {
                    // if cut inequality is violated

                    // do a DFS to find all the edges, that leave the partition, which is essentially a cut
                    auto label = components_.labels_[lv0];

                    std::fill(visited_.begin(), visited_.end(), 0);

                    std::stack<size_t> S;
                    S.push(lv0);
                    visited_[lv0] = 1;

                    ptrdiff_t sz = 0;
                    while (!S.empty())
                    {
                        auto const v = S.top();
                        S.pop();

                        for (auto it = original_graph_.adjacenciesFromVertexBegin(v); it != original_graph_.adjacenciesFromVertexEnd(v); ++it)
                            if (components_.labels_[it->vertex()] != label)
                            {
                                coefficients_[sz] = -1;
                                variables_[sz] = edge_in_lifted_graph_[it->edge()];
                                
                                ++sz;
                            }
                            else if (!visited_[it->vertex()])
                            {
                                S.push(it->vertex());
                                visited_[it->vertex()] = 1;
                            }
                    }

                    coefficients_[sz] = 1;
                    variables_[sz] = edge;

                    this->addLazyConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                    ++nCut;


                    label = components_.labels_[lv1];

                    S.push(lv1);
                    visited_[lv1] = 1;

                    sz = 0;
                    while (!S.empty())
                    {
                        auto const v = S.top();
                        S.pop();

                        for (auto it = original_graph_.adjacenciesFromVertexBegin(v); it != original_graph_.adjacenciesFromVertexEnd(v); ++it)
                            if (components_.labels_[it->vertex()] != label)
                            {
                                coefficients_[sz] = -1;
                                variables_[sz] = edge_in_lifted_graph_[it->edge()];
                                
                                ++sz;
                            }
                            else if (!visited_[it->vertex()])
                            {
                                S.push(it->vertex());
                                visited_[it->vertex()] = 1;
                            }
                    }

                    coefficients_[sz] = 1;
                    variables_[sz] = edge;

                    this->addLazyConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                    ++nCut;
                }
            }

            std::cerr << std::fixed << std::setprecision(10) << this->objectiveBound_ << " " << std::setprecision(10) << this->objectiveBest_ << " " << nCycle << " " << nPath << " " << nCut << std::endl;
        }

    private:
        struct SubgraphWithCut
        {
            SubgraphWithCut(Callback& callback, std::vector<size_t> const& edge_in_lifted_graph) :
                callback_(callback), edge_in_lifted_graph_(edge_in_lifted_graph)
            {}

            bool vertex(size_t v) const
            {
                return true;
            }

            bool edge(size_t e) const
            {
                return callback_.label(edge_in_lifted_graph_[e]) < .5;
            }

            Callback& callback_;
            std::vector<size_t> const& edge_in_lifted_graph_;
        };

        ORIGGRAPH const& original_graph_;
        LIFTGRAPH const& lifted_graph_;

        std::vector<double> coefficients_;
        ComponentsBySearch<ORIGGRAPH> components_;
        std::vector<ptrdiff_t> buffer_;
        std::vector<size_t> edge_in_lifted_graph_;
        std::deque<size_t> path_;
        std::vector<size_t> variables_;
        std::vector<char> visited_;
    };

    ILP ilp;

    ilp.setRelativeGap(0.0);
    ilp.setAbsoluteGap(0.0);
    ilp.setTimeLimit(timeLimitSeconds);
    ilp.addVariables(edgeCosts.size(), edgeCosts.data());

    Callback callback(ilp, original_graph, lifted_graph);
    ilp.setCallback(callback);

    ilp.optimize();

    std::vector<size_t> edge_in_lifted_graph(original_graph.numberOfEdges());
    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto const v0 = original_graph.vertexOfEdge(i, 0);
        auto const v1 = original_graph.vertexOfEdge(i, 1);

        edge_in_lifted_graph[i] = lifted_graph.findEdge(v0, v1).second;
    }

    ComponentsBySearch<ORIGGRAPH> components;
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

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_ILP_CALLBACK_HXX
