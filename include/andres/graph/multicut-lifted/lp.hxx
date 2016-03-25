#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_LP_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_LP_HXX

#include <vector>
#include <deque>
#include <algorithm>
#include <stack>

#include <andres/graph/shortest-paths.hxx>
#include <levinkov/timer.hxx>


namespace andres {
namespace graph {
namespace multicut_lifted {

template<typename LP, typename ORIGGRAPH, typename LIFTGRAPH, typename ECA>
inline
std::vector<double> lp(ORIGGRAPH const& original_graph, LIFTGRAPH const& lifted_graph, ECA const& edgeCosts, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct Visitor
    {
        bool operator()() const
        {
            return true;
        }
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
        DinicFlow(int n) :
            dist_(n), q_(n), work_(n), g_(n)
        {}

        // Adds bidirectional edge
        void addEdge(int s, int t, int cap)
        {
            Edge a = { t, static_cast<int>(g_[t].size()), 0, cap };
            Edge b = { s, static_cast<int>(g_[s].size()), 0, cap };

            g_[s].push_back(a);
            g_[t].push_back(b);
        }

        void clear()
        {
            for (auto& v : g_)
                v.clear();
        }

        int maxFlow(int src, int dest)
        {
            src_ = src;
            dest_ = dest;
            
            int result = 0;

            while (bfs())
            {
                std::fill(work_.begin(), work_.end(), 0);

                while (auto delta = dfs(src_, std::numeric_limits<int>::max()))
                    result += delta;
            }

            return result;
        }

        // return edges (as pairs of vertices) of the Min Cut
        std::set<std::pair<int, int>> getMinCut()
        {
            std::fill(work_.begin(), work_.end(), 0);

            std::stack<int> S;

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

            std::set<std::pair<int, int>> ret;

            for (int i = 0; i < g_.size(); ++i)
                for (auto& e : g_[i])
                    if (work_[i] != work_[e.to])
                        ret.insert(std::make_pair(std::min(i, e.to), std::max(i, e.to)));

            return ret;
        }

    private:
        struct Edge
        {
            int to, rev;
            int f, cap;
        };

        bool bfs()
        {
            std::fill(dist_.begin(), dist_.end(), -1);

            dist_[src_] = 0;

            int qt = 0;
            q_[qt++] = src_;

            for (int qh = 0; qh < qt; qh++)
            {
                auto u = q_[qh];

                for (int j = 0; j < g_[u].size(); j++)
                {
                    Edge &e = g_[u][j];
                    auto v = e.to;
                    
                    if (dist_[v] < 0 && e.f < e.cap)
                    {
                        dist_[v] = dist_[u] + 1;
                        q_[qt++] = v;
                    }
                }
            }

            return dist_[dest_] >= 0;
        }

        int dfs(int u, int f)
        {
            if (u == dest_)
                return f;

            for (int &i = work_[u]; i < (int) g_[u].size(); i++)
            {
                Edge &e = g_[u][i];

                if (e.cap <= e.f) 
                    continue;

                int v = e.to;

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

        int src_, dest_;
        std::vector<int> dist_, q_, work_;
        std::vector<std::vector<Edge>> g_;
    };

    LP lp;
    levinkov::Timer t;
    
    std::vector<double> coefficients(lifted_graph.numberOfEdges());
    std::vector<double> distances(lifted_graph.numberOfVertices()); 
    std::vector<size_t> edge_in_lifted_graph(original_graph.numberOfEdges());
    std::vector<size_t> parents(lifted_graph.numberOfVertices());
    std::deque<size_t> path;
    std::vector<double> vars(original_graph.numberOfEdges());
    std::vector<size_t> variables(lifted_graph.numberOfEdges());

    auto separateAndAddInequalities = [&] ()
    {
        t.stop();

        levinkov::Timer t_separation;
        t_separation.start();

        DinicFlow flow(original_graph.numberOfVertices());

        for (size_t i = 0; i < vars.size(); ++i)
        {
            vars[i] = std::min(std::max(.0, lp.variableValue(edge_in_lifted_graph[i])), 1.0);

            auto v0 = lifted_graph.vertexOfEdge(edge_in_lifted_graph[i], 0);
            auto v1 = lifted_graph.vertexOfEdge(edge_in_lifted_graph[i], 1);

            flow.addEdge(v0, v1, 100000.0*(1.0 - vars[i]));
        }

        // search for violated non-chordal cycles and add corresp. inequalities
        size_t nCycle = 0;
        size_t nPath = 0;
        size_t nCut = 0;

        for (size_t edge = 0; edge < lifted_graph.numberOfEdges(); ++edge)
        {
            auto lv0 = lifted_graph.vertexOfEdge(edge, 0);
            auto lv1 = lifted_graph.vertexOfEdge(edge, 1);

            // search for shortest path
            double distance;
            spsp(original_graph, DefaultSubgraphMask<>(), lv0, lv1, vars.begin(), path, distance, distances.begin(), parents.begin());
            
            if (lp.variableValue(edge) > distance + tolerance)
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

            flow.maxFlow(lv0, lv1);

            double S_cut = .0;

            ptrdiff_t sz = 0;
            for (auto& p : flow.getMinCut())
            {
                coefficients[sz] = -1.0;
                variables[sz] = lifted_graph.findEdge(p.first, p.second).second;

                S_cut += 1.0 - vars[variables[sz]];

                ++sz;
            }

            if (1.0 - lp.variableValue(edge) > S_cut + tolerance)
            {
                coefficients[sz] = 1.0;
                variables[sz] = edge;

                lp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());    

                ++nCut;
            }
        }

        t_separation.stop();

        double objValue = .0;
        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
            objValue += lp.variableValue(i)*edgeCosts[i];

        std::cerr << t.get_elapsed_seconds() << " " << objValue << " " << nCycle << " " << nPath << " " << nCut << " " << t_separation.get_elapsed_seconds() << std::endl;

        t.start();

        return nCycle + nPath + nCut;
    };

    t.start();

    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto v0 = original_graph.vertexOfEdge(i, 0);
        auto v1 = original_graph.vertexOfEdge(i, 1);

        edge_in_lifted_graph[i] = lifted_graph.findEdge(v0, v1).second;
    }

    lp.initModel(lifted_graph.numberOfEdges(), edgeCosts.data());

    for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        if (!visitor())
            break;

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
