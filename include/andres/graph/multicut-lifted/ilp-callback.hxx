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

// Solver for the minimum cost lifted multicut problem for arbitrary graphs.
//
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
        SubgraphWithCut(ILP const& ilp, std::vector<size_t> const& edge_index_lifted) :
            ilp_(ilp), edge_index_lifted_(edge_index_lifted)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return ilp_.label(edge_index_lifted_[e]) < .5;
        }

        ILP const& ilp_;
        std::vector<size_t> const& edge_index_lifted_;
    };

    // Dinic's Max Flow algorithm for undirected graphs
    class DinicFlow
    {
    public:
        DinicFlow(size_t n) :
            dist_(n), q_(n), work_(n), g_(n)
        {}

        // Adds bidirectional edge
        void addEdge(size_t s, size_t t, ptrdiff_t cap)
        {
            Edge a = { t, g_[t].size(), 0, cap };
            Edge b = { s, g_[s].size(), 0, cap };

            g_[s].push_back(a);
            g_[t].push_back(b);
        }

        void clear()
        {
            for (auto& v : g_)
                v.clear();
        }

        size_t maxFlow(size_t src, size_t dest)
        {
            src_ = src;
            dest_ = dest;

            for (auto& adj : g_)
                for (auto& e : adj)
                    e.f = 0;

            size_t result = 0;

            while (bfs())
            {
                std::fill(work_.begin(), work_.end(), 0);

                while (auto delta = dfs(src_, std::numeric_limits<ptrdiff_t>::max()))
                    result += delta;
            }

            return result;
        }

        // return edges (as pairs of vertices) of the Min Cut
        std::set<std::pair<size_t, size_t>> getMinCut()
        {
            std::fill(work_.begin(), work_.end(), 0);

            std::stack<size_t> S;

            work_[src_] = 1;
            S.push(src_);

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                for (auto& e : g_[v])
                    if (e.f < e.cap && work_[e.to] == 0)
                    {
                        work_[e.to] = 1;
                        S.push(e.to);
                    }
            }

            std::set<std::pair<size_t, size_t>> ret;

            for (size_t i = 0; i < g_.size(); ++i)
                for (auto& e : g_[i])
                    if (work_[i] != work_[e.to])
                        ret.insert(std::make_pair(std::min(i, e.to), std::max(i, e.to)));

            return ret;
        }

    private:
        struct Edge
        {
            size_t to, rev;
            ptrdiff_t f, cap;
        };

        bool bfs()
        {
            std::fill(dist_.begin(), dist_.end(), -1);

            dist_[src_] = 0;

            size_t qt = 0;
            q_[qt++] = src_;

            for (size_t qh = 0; qh < qt; qh++)
            {
                auto u = q_[qh];

                for (size_t j = 0; j < g_[u].size(); j++)
                {
                    auto& e = g_[u][j];
                    auto  v = e.to;
                    
                    if (dist_[v] < 0 && e.f < e.cap)
                    {
                        dist_[v] = dist_[u] + 1;
                        q_[qt++] = v;
                    }
                }
            }

            return dist_[dest_] >= 0;
        }

        size_t dfs(size_t u, ptrdiff_t f)
        {
            if (u == dest_)
                return f;

            for (auto& i = work_[u]; i < g_[u].size(); i++)
            {
                auto& e = g_[u][i];

                if (e.cap <= e.f) 
                    continue;

                auto v = e.to;

                if (dist_[v] == dist_[u] + 1)
                {
                    auto df = dfs(v, std::min(f, e.cap - e.f));

                    if (df > 0)
                    {
                        e.f += df;
                        g_[v][e.rev].f -= df;

                        return df;
                    }
                }
            }

            return 0;
        }

        size_t src_, dest_;
        std::vector<int> dist_;
        std::vector<size_t> q_, work_;
        std::vector<std::vector<Edge>> g_;
    };
    
    class Callback: public ILP::Callback
    {
    public:
        Callback(ILP& solver, ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts) :
            ILP::Callback(solver),
            original_graph_(original_graph),
            lifted_graph_(lifted_graph),
            coefficients_(lifted_graph.numberOfEdges()),
            variables_(lifted_graph.numberOfEdges()),
            visited_(original_graph.numberOfVertices()),
            edge_index_lifted_(original_graph.numberOfEdges()),
            edgeCosts_(edgeCosts)
        {
            // find original edges in lifted graph
            for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
            {
                auto const v0 = original_graph.vertexOfEdge(i, 0);
                auto const v1 = original_graph.vertexOfEdge(i, 1);

                edge_index_lifted_[i] = lifted_graph.findEdge(v0, v1).second;
            }
        }

        void separateAndAddLazyConstraints() override
        {
            // connected component labeling of the original graph
            auto nComponents = components_.build(original_graph_, SubgraphWithCut(*this, edge_index_lifted_));

            // MinCut setup
            DinicFlow flow(original_graph_.numberOfVertices());
            std::vector<double> vars(original_graph_.numberOfEdges());
            for (size_t i = 0; i < vars.size(); ++i)
            {
                vars[i] = this->label(edge_index_lifted_[i]);

                auto const v0 = lifted_graph_.vertexOfEdge(edge_index_lifted_[i], 0);
                auto const v1 = lifted_graph_.vertexOfEdge(edge_index_lifted_[i], 1);

                if (vars[i] > .5)
                    flow.addEdge(v0, v1, 1);
                else
                    flow.addEdge(v0, v1, original_graph_.numberOfEdges());
            }

            // search for violated non-chordal cycles and add corresp. inequalities
            size_t nCycle = 0;
            size_t nPath = 0;
            size_t nCut = 0;
            size_t nGenCycle = 0;

            std::vector<ptrdiff_t> parent(lifted_graph_.numberOfVertices());
            std::vector<double> cost(lifted_graph_.numberOfVertices());
            std::vector<size_t> dist(lifted_graph_.numberOfVertices());

            // iterate over all edges of the super graph
            for (ptrdiff_t edge = 0; edge < lifted_graph_.numberOfEdges(); ++edge)
            {
                auto const lv0 = lifted_graph_.vertexOfEdge(edge, 0);
                auto const lv1 = lifted_graph_.vertexOfEdge(edge, 1);

                // edge violates cycle/path inequality
                if (this->label(edge) > .5 && components_.areConnected(lv0, lv1))
                {
                    // The following is a specialized implementation of a bidirectional BFS to compute a shortest
                    // path from lv0 to lv1. The search is based on the connectivity defined by the connected
                    // component labeling of the original graph but takes shortcuts via edges of the super graph.

                    // The current version also maximizes the sum of edge costs among all shortest (# of hops) paths
                    // (this roughly doubles the number of explored nodes in the BFS, but increases the strength
                    // of the separated inequalities wrt to the costs of the objective).

                    std::fill(parent.begin(), parent.end(), 0);
                    std::fill(cost.begin(), cost.end(), 0);
                    std::fill(dist.begin(), dist.end(), 0);

                    std::queue<size_t> queues[2];
                    queues[0].push(lv0);
                    queues[1].push(lv1);

                    parent[lv0] = lv0 + 1;
                    parent[lv1] = -static_cast<std::ptrdiff_t>(lv1) - 1;

                    size_t curr_dist = 0;

                    // vertices corresponding to optimal search tree connection
                    size_t v0;
                    size_t v1;
                    
                    bool found_path = false;
                    double path_cost = -std::numeric_limits<double>::infinity();
                        
                    for (size_t q = 0; !found_path && (!queues[0].empty() || !queues[1].empty()); q = 1 - q)
                    {
                        if (q == 0)
                            curr_dist++;

                        const size_t nFrontierNodes = queues[q].size();
                        for (size_t n = 0; n < nFrontierNodes; n++)
                        {
                            auto v = queues[q].front();
                            queues[q].pop();

                            for (auto it = lifted_graph_.adjacenciesFromVertexBegin(v); it != lifted_graph_.adjacenciesFromVertexEnd(v); it++)
                            {
                                auto w = it->vertex();
                                auto e = it->edge();

                                // find unexplored vertex
                                if (this->label(e) < .5 && components_.areConnected(lv0,w) && !parent[w])
                                {
                                    queues[q].push(w);

                                    if (q == 0)
                                        parent[w] = v + 1;
                                    else
                                        parent[w] = -static_cast<std::ptrdiff_t>(v) - 1;

                                    dist[w] = curr_dist;
                                    cost[w] = cost[v] + edgeCosts_[e];

                                }
                                // update path cost to frontier node
                                else if (this->label(e) < .5 && q == 0 && parent[w] > 0 && dist[w] == curr_dist && cost[w] < cost[v] + edgeCosts_[e])
                                {
                                    parent[w] = v + 1;
                                    cost[w] = cost[v] + edgeCosts_[e];
                                }
                                else if (this->label(e) < .5 && q == 1 && parent[w] < 0 && dist[w] == curr_dist && cost[w] < cost[v] + edgeCosts_[e])
                                {
                                    parent[w] = -static_cast<std::ptrdiff_t>(v) - 1;
                                    cost[w] = cost[v] + edgeCosts_[e];
                                }
                                // if the two search trees connect
                                else if (this->label(e) < .5 && parent[w] != 0 && cost[v] + cost[w] + edgeCosts_[e] > path_cost)
                                {
                                    if (q == 0 && parent[w] < 0)
                                    {
                                        v0 = v;
                                        v1 = w;
                                    }
                                    else if (q == 1 && parent[w] > 0)
                                    {
                                        v0 = w;
                                        v1 = v;
                                    }
                                    else
                                        continue;

                                    path_cost = cost[v] + cost[w] + edgeCosts_[e];
                                    found_path = true;
                                }                                
                            }
                        }
                    }

                    // DEBUG:
                    if (!found_path)
                    {
                        std::cout << "Warning: no path found for edge " << edge << std::endl;
                        continue; 
                    }               

                    // fill path from bidirectional BFS
                    path_.clear();
                    auto p = v0;
                    path_.push_front(p);
                    while (parent[p] - 1 != p)
                    {
                        p = parent[p] - 1;
                        path_.push_front(p);
                    }
                    p = v1;
                    path_.push_back(p);
                    while (-parent[p] - 1 != p)
                    {
                        p = -parent[p] - 1;
                        path_.push_back(p);
                    }

                    // skip chord check for triangles
                    if (path_.size() > 3)
                    {
                        // check for chords in linear time
                        std::fill(visited_.begin(), visited_.end(), 0);
                        if (hasChord(lifted_graph_, path_.begin(), path_.end(), visited_, true))
                            continue;
                    }

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
                // edge violates cut inequality
                else if (this->label(edge) < .5 && !components_.areConnected(lv0, lv1))
                {
                    // search for triangle inequality with lifted edge on right hand side and best associated cost
                    bool add_triangle = false;
                    double best_cost = std::numeric_limits<double>::infinity();

                    for (auto it1 = lifted_graph_.adjacenciesFromVertexBegin(lv0); it1 != lifted_graph_.adjacenciesFromVertexEnd(lv0); it1++)
                    {
                        auto v1 = it1->vertex();

                        auto e2p = lifted_graph_.findEdge(lv1,v1);

                        if (v1 != lv0 && e2p.first)
                        {
                            auto e1 = it1->edge();
                            auto e2 = e2p.second;

                            if (this->label(e1) < .5 && this->label(e2) > .5 && edgeCosts_[e2] - edgeCosts_[e1] < best_cost)
                            {
                                best_cost = edgeCosts_[e1] - edgeCosts_[e2];

                                variables_[0] = edge;
                                variables_[1] = e1;
                                variables_[2] = e2;

                                add_triangle = true;
                                break;
                            }
                            else if (this->label(e1) > .5 && this->label(e2) < .5 && edgeCosts_[e1] - edgeCosts_[e2] < best_cost)
                            {
                                best_cost = edgeCosts_[e2] - edgeCosts_[e1];

                                variables_[0] = edge;
                                variables_[1] = e2;
                                variables_[2] = e1;

                                add_triangle = true;
                                break;
                            }
                        }
                        
                    }

                    // if triangle inequality is found, skip search for cut inequality
                    if (add_triangle)
                    {
                        coefficients_[0] = 1.0;
                        coefficients_[1] = 1.0;
                        coefficients_[2] = -1.0;

                        this->addLazyConstraint(variables_.begin(), variables_.begin() + 3, coefficients_.begin(), 0, std::numeric_limits<double>::infinity());    
                        ++nGenCycle;
                        continue;
                    }                        

                    // search for minimum size cut-set separating lv0 and lv1
                    ptrdiff_t sz = 0;
                    auto flow_value = static_cast<double>(flow.maxFlow(lv0, lv1));
                    
                    for (auto& p : flow.getMinCut())
                    {
                        coefficients_[sz] = -1.0;                            
                        
                        variables_[sz] = lifted_graph_.findEdge(p.first, p.second).second;
                        ++sz;
                    }
 
                    coefficients_[sz] = 1.0;
                    variables_[sz] = edge;

                    this->addLazyConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                    ++nCut;
                }
            }

            std::cout << this->getDoubleInfo(GRB_CB_RUNTIME) << " " << this->getDoubleInfo(GRB_CB_MIPSOL_NODCNT) << " " << std::flush;
            std::cerr << std::fixed << std::setprecision(5) << this->objectiveBound_ << " " << this->objectiveBest_ << " " << nCycle << " " << nPath << " " << nCut << " " << nGenCycle << std::endl;
        }

    private:
        struct SubgraphWithCut
        {
            SubgraphWithCut(Callback& callback, std::vector<size_t> const& edge_index_lifted) :
                callback_(callback), edge_index_lifted_(edge_index_lifted)
            {}

            bool vertex(size_t v) const
            {
                return true;
            }

            bool edge(size_t e) const
            {
                return callback_.label(edge_index_lifted_[e]) < .5;
            }

            Callback& callback_;
            std::vector<size_t> const& edge_index_lifted_;
        };

        ORIGGRAPH const& original_graph_;
        LIFTGRAPH const& lifted_graph_;
        ECA const& edgeCosts_;

        std::vector<double> coefficients_;
        ComponentsBySearch<ORIGGRAPH> components_;
        std::vector<ptrdiff_t> buffer_;
        std::vector<size_t> edge_index_lifted_;
        std::deque<size_t> path_;
        std::vector<size_t> variables_;
        std::vector<char> visited_;
    };

    ILP ilp;

    ilp.setRelativeGap(0.0);
    ilp.setAbsoluteGap(0.0);
    ilp.setTimeLimit(timeLimitSeconds);
    ilp.addVariables(edgeCosts.size(), edgeCosts.data());

    Callback callback(ilp, original_graph, lifted_graph, edgeCosts);
    ilp.setCallback(callback);

    ilp.optimize();

    std::vector<size_t> edge_index_lifted(original_graph.numberOfEdges());
    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto const v0 = original_graph.vertexOfEdge(i, 0);
        auto const v1 = original_graph.vertexOfEdge(i, 1);

        edge_index_lifted[i] = lifted_graph.findEdge(v0, v1).second;
    }

    ComponentsBySearch<ORIGGRAPH> components;
    components.build(original_graph, SubgraphWithCut(ilp, edge_index_lifted));

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
