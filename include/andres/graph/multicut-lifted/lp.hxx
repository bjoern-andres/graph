#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_LP_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_LP_HXX

#include <vector>
#include <deque>
#include <algorithm>
#include <stack>

#include <andres/graph/shortest-paths.hxx>


namespace andres {
namespace graph {
namespace multicut_lifted {

template<typename LP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA>
inline
std::vector<double> lp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct EmptyVisitor
    {
    } visitor;

    return lp<LP>(original_graph, lifted_graph, edgeCosts, visitor, numberOfIterations);
}

template<typename LP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename VIS>
inline
std::vector<double> lp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, VIS& visitor, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    const double tolerance = std::numeric_limits<float>::epsilon();

    // Dinic's Max FLow algorithm for undirected graph
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

    LP lp;
    
    std::vector<double> coefficients(lifted_graph.numberOfEdges());
    std::vector<double> distances(lifted_graph.numberOfVertices()); 
    std::vector<size_t> edge_in_lifted_graph(original_graph.numberOfEdges());
    std::vector<size_t> parents(lifted_graph.numberOfVertices());
    std::deque<size_t> path;
    std::vector<double> vars(original_graph.numberOfEdges());
    std::vector<size_t> variables(lifted_graph.numberOfEdges());

    auto separateAndAddInequalities = [&] ()
    {
        DinicFlow flow(original_graph.numberOfVertices());

        for (size_t i = 0; i < vars.size(); ++i)
        {
            vars[i] = std::min(std::max(.0, lp.variableValue(edge_in_lifted_graph[i])), 1.0);

            auto const v0 = lifted_graph.vertexOfEdge(edge_in_lifted_graph[i], 0);
            auto const v1 = lifted_graph.vertexOfEdge(edge_in_lifted_graph[i], 1);

            // as MaxLow can be computed only for intergral edge weights, we lose a bit of precision doing the following
            // if your weights are very small, you might want to increase the constant (everywhere in this file),
            // in order not to run into precision problems
            flow.addEdge(v0, v1, 100000000000.0*(1.0 - vars[i]));
        }

        // search for violated non-chordal cycles and add corresp. inequalities
        size_t nCycle = 0;
        size_t nPath = 0;
        size_t nCut = 0;

        for (size_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
        {
            auto const lv0 = lifted_graph.vertexOfEdge(edge, 0);
            auto const lv1 = lifted_graph.vertexOfEdge(edge, 1);

            // search for shortest path
            double distance;
            spsp(original_graph, DefaultSubgraphMask<>(), lv0, lv1, vars.begin(), path, distance, distances.begin(), parents.begin());

            bool chordless = true;
            for (auto it1 = path.begin(); it1 != path.end() - 2 && chordless; ++it1)
                for (auto it2 = it1 + 2; it2 != path.end(); ++it2)
                {
                    if (it1 == path.begin() && it2 == path.end() - 1)
                        continue;

                    auto const e = lifted_graph.findEdge(*it1, *it2);
                    if (e.first && std::min(std::max(.0, lp.variableValue(e.second)), 1.0) > distances[*it2] - distances[*it1] + tolerance)
                    {
                        chordless = false;
                        break;
                    }
                }

            if (!chordless)
                continue;
            
            if (std::min(std::max(.0, lp.variableValue(edge)), 1.0) > distance + tolerance)
            {
                // add inequality
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    coefficients[j] = 1.0;
                    variables[j] = lifted_graph.findEdge(path[j], path[j + 1]).second;
                }

                coefficients[path.size() - 1] = -1.0;
                variables[path.size() - 1] = edge;

                lp.addConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());

                if (original_graph.findEdge(lv0, lv1).first)
                    ++nCycle;
                else
                    ++nPath;
            }

            if (!original_graph.findEdge(lv0, lv1).first)
            {
                // find min cut only for lifted edges
                // find cut(s) for both vertices of an edge - considerably improves convergence

                auto flow_value = static_cast<double>(flow.maxFlow(lv0, lv1)) / 100000000000.0;

                if (1.0 - std::min(std::max(.0, lp.variableValue(edge)), 1.0) > flow_value + tolerance)
                {
                    ptrdiff_t sz = 0;
                    for (auto& p : flow.getMinCut())
                    {
                        coefficients[sz] = -1.0;
                        variables[sz] = lifted_graph.findEdge(p.first, p.second).second;

                        ++sz;
                    }

                    coefficients[sz] = 1.0;
                    variables[sz] = edge;

                    lp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                    ++nCut;
                }

                flow_value = static_cast<double>(flow.maxFlow(lv1, lv0)) / 100000000000.0;

                if (1.0 - std::min(std::max(.0, lp.variableValue(edge)), 1.0) > flow_value + tolerance)
                {
                    ptrdiff_t sz = 0;
                    for (auto& p : flow.getMinCut())
                    {
                        coefficients[sz] = -1.0;
                        variables[sz] = lifted_graph.findEdge(p.first, p.second).second;

                        ++sz;
                    }

                    coefficients[sz] = 1.0;
                    variables[sz] = edge;

                    lp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                    ++nCut;
                }
            }
        }

        double objValue = .0;
        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
            objValue += lp.variableValue(i)*edgeCosts[i];

        std::cerr << std::fixed << std::setprecision(10) << objValue << " " << nCycle << " " << nPath << " " << nCut << std::endl;

        return nCycle + nPath + nCut;
    };

    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto const v0 = original_graph.vertexOfEdge(i, 0);
        auto const v1 = original_graph.vertexOfEdge(i, 1);

        edge_in_lifted_graph[i] = lifted_graph.findEdge(v0, v1).second;
    }

    lp.initModel(lifted_graph.numberOfEdges(), edgeCosts.data());

    for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        lp.optimize();

        if (separateAndAddInequalities() == 0)
            break;
    }

    std::vector<double> edge_values(lifted_graph.numberOfEdges());
    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        edge_values[i] = lp.variableValue(i);

    return edge_values;
}

}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_LP_HXX
