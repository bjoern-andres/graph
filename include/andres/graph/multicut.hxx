#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_HXX
#define ANDRES_GRAPH_MULTICUT_HXX

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <deque>
#include <algorithm> // std::copy

#include "andres/partition.hxx" 
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/paths.hxx"
#include "andres/graph/components.hxx"
#include "andres/graph/shortest-paths.hxx"

namespace andres {
namespace graph {

/// Solver for the Minimum Cost Multicut Problem for arbitrary graphs.
///
/// This is a variant of the solver proposed in 
/// 
/// Andres B., Kroeger T., Briggman K. L., Denk W., Korogod N., Knott G., Koethe U. and Hamprecht F. A.
/// Globally Optimal Closed-surface Segmentation for Connectomics. ECCV 2012
/// http://dx.doi.org/10.1007/978-3-642-33712-3_56
///
/// This code operates on graphs whereas the code used to produce the results 
/// in the above publication operates on cellular complexes. While cellular 
/// complexes are richer structures than graphs and facilitates additional 
/// tweaks of the solver (e.g. cell suppression), graphs are more common. 
/// This code has a wider range of applications.
/// 
template<class GRAPH, class ILP>
class Multicut {
public:
    typedef GRAPH GraphType;
    typedef ILP Ilp;
    typedef ComponentsBySearch<GRAPH> Components;
    
    /// Visitors can be used to follow the progress of the optimization.
    struct IdleVisitor {
        void repairSolution() const {}
        void optimize() const {}
        void addCycleInequalities() const {}
        void endIteration(const std::size_t, const std::size_t) const {}
    };

    Multicut();
    void setup(const GraphType&, const std::vector<double>&);
    void solve(const std::size_t);
    template<class VISITOR>
        void solve(const std::size_t, VISITOR&);
    
    Ilp& ilp();
    double label(const std::size_t) const;

private:
    struct SubgraphWithoutCut { // a subgraph mask
        SubgraphWithoutCut(const ILP& ilp)
            : ilp_(ilp) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return ilp_.label(e) == 0; }

        const ILP& ilp_;
    };

    std::size_t addCycleInequalities();
    void repairSolution();

    Components components_;
    const GraphType* graph_;
    Ilp ilp_;
};

/// Solver for the Minimum Cost Multicut Problem for complete graphs (Set Partition Problem).
template<class GRAPH_VISITOR, class ILP>
class Multicut<CompleteGraph<GRAPH_VISITOR>, ILP> {
public:
    typedef CompleteGraph<GRAPH_VISITOR> GraphType;
    typedef ILP Ilp;
    typedef ComponentsBySearch<CompleteGraph<GRAPH_VISITOR> > Components;

    /// Visitors can be used to follow the progress of the optimization.
    struct IdleVisitor {
        void repairSolution() const {}
        void optimize() const {}
        void addCycleInequalities() const {}
        void endIteration(const std::size_t, const std::size_t) const {}
    };

    Multicut();
    void setup(const GraphType&, const std::vector<double>&);
    void solve(const std::size_t);
    template<class VISITOR>
        void solve(const std::size_t, VISITOR&);

    Ilp& ilp();
    double label(const std::size_t) const;

private:
    struct SubgraphWithoutCut { // a subgraph mask
        SubgraphWithoutCut(const ILP& ilp)
            : ilp_(ilp) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return ilp_.label(e) == 0; }

        const ILP& ilp_;
    };

    std::size_t addCycleInequalities();
    void repairSolution();

    Components components_;
    const GraphType* graph_;
    Ilp ilp_;
};

template<class GRAPH, class ILP>
inline
Multicut<GRAPH, ILP>::Multicut()
:   components_(),
    graph_(NULL),
    ilp_()
{}

template<class GRAPH, class ILP>
inline double
Multicut<GRAPH, ILP>::label(
    const std::size_t edge
) const {
    assert(edge < graph_->numberOfEdges());
    size_t v0 = graph_->vertexOfEdge(edge, 0);
    size_t v1 = graph_->vertexOfEdge(edge, 1);
    return components_.areConnected(v0, v1) ? 0 : 1;
}

template<class GRAPH, class ILP>
inline typename Multicut<GRAPH, ILP>::Ilp&
Multicut<GRAPH, ILP>::ilp() {
    return ilp_;
}

template<class GRAPH, class ILP>
void 
Multicut<GRAPH, ILP>::setup(
    const GraphType& graph,
    const std::vector<double>& edgeCosts
) {
    if(graph.numberOfEdges() != edgeCosts.size()) {
        throw std::runtime_error("the number of edges in the graph does not equal the number of costs.");
    }
    graph_ = &graph;
    ilp_.initModel(graph_->numberOfEdges(), edgeCosts.data());
}

template<class GRAPH, class ILP>
inline void 
Multicut<GRAPH, ILP>::solve(
    const std::size_t maxIterations
) {
    IdleVisitor idleVisitor;
    solve(maxIterations, idleVisitor);
}

template<class GRAPH, class ILP>
template<class VISITOR>
void 
Multicut<GRAPH, ILP>::solve(
    const std::size_t maxIterations,
    VISITOR& visitor
) {
    std::size_t i = 0;

    for (; maxIterations == 0 || i < maxIterations; ++i) {
        if(i != 0) {
            visitor.repairSolution();
            repairSolution();
        }
        visitor.optimize();
        ilp_.optimize();
        visitor.addCycleInequalities();
        std::size_t nc = addCycleInequalities();
        visitor.endIteration(i, nc);        
        if(nc == 0) {
            break;
        }
    }

    if (i == maxIterations)
    {
        visitor.repairSolution();
        repairSolution();
    }
}

template<class GRAPH, class ILP>
inline std::size_t
Multicut<GRAPH, ILP>::addCycleInequalities() {
    // label connected components
    components_.build(*graph_, SubgraphWithoutCut(ilp_));

    // search for violated non-chordal cycles and add corresp. inequalities
    std::deque<std::size_t> path;
    std::vector<ptrdiff_t> buffer;
    std::vector<double> variables;
    std::vector<double> coefficients;
    std::size_t counter = 0;
    #pragma omp parallel for private(path, buffer, variables, coefficients), schedule(guided)
    for(ptrdiff_t edge = 0; edge < graph_->numberOfEdges(); ++edge) {
        if(ilp_.label(edge) == 1) {
            const std::size_t v0 = graph_->vertexOfEdge(edge, 0);
            const std::size_t v1 = graph_->vertexOfEdge(edge, 1);
            if(components_.areConnected(v0, v1)) { 
                // search for shortest path
                bool found = spsp(*graph_, SubgraphWithoutCut(ilp_), v0, v1, path, buffer); 
                assert(found);
                
                // skip chordal paths
                std::pair<bool, std::size_t> chordal;
                chordal = findChord(*graph_, path.begin(), path.end(), true);
                if(chordal.first) {
                    continue;
                }

                // add inequality
                variables.resize(path.size());
                for(std::size_t j = 0; j < path.size() - 1; ++j) {
                    std::pair<bool, std::size_t> p = graph_->findEdge(path[j], path[j + 1]);
                    assert(p.first);
                    variables[j] = static_cast<double>(p.second);
                }              
                variables.back() = static_cast<double>(edge);
                coefficients.resize(variables.size());
                std::fill(coefficients.begin(), coefficients.end(), 1);
                coefficients.back() = -1;
                #pragma omp critical
                ilp_.addConstraint(
                    variables.begin(), variables.end(),
                    coefficients.begin(),
                    0, std::numeric_limits<double>::infinity()
                );

                #pragma omp atomic
                ++counter;
            }
        }
    }
    return counter;
}

// pre-condition: componenets_ must be up to date
template<class GRAPH, class ILP>
inline void
Multicut<GRAPH, ILP>::repairSolution() {
    std::vector<double> repairedSolution(graph_->numberOfEdges());
    for(std::size_t edge = 0; edge < graph_->numberOfEdges(); ++edge) {
        const std::size_t v0 = graph_->vertexOfEdge(edge, 0);
        const std::size_t v1 = graph_->vertexOfEdge(edge, 1);
        if(!components_.areConnected(v0, v1)) {
            repairedSolution[edge] = 1;
        }
    }
    ilp_.setStart(repairedSolution.begin());
}

template<class GRAPH_VISITOR, class ILP>
inline
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::Multicut()
:   components_(),
    graph_(NULL),
    ilp_()
{}

template<class GRAPH_VISITOR, class ILP>
inline void
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::setup(
    const GraphType& graph,
    const std::vector<double>& edgeCosts
) {
    if(graph.numberOfEdges() != edgeCosts.size()) {
        throw std::runtime_error("the number of edges in the graph does not equal the number of costs.");
    }
    graph_ = &graph;
    ilp_.initModel(graph_->numberOfEdges(), edgeCosts.data());
}

template<class GRAPH_VISITOR, class ILP>
inline void
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::solve(
    const std::size_t maxIterations
) {
    IdleVisitor idleVisitor;
    solve(maxIterations, idleVisitor);
}

template<class GRAPH_VISITOR, class ILP>
template<class VISITOR>
void
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::solve(
    const std::size_t maxIterations,
    VISITOR& visitor
) {
    std::size_t i = 0;

    for(; maxIterations == 0 || i < maxIterations; ++i) {
        // not repairing solution because connected component labeling
        // (whose runtime complexity is linear in the number of nodes plus the
        // number of edges) is usually slower for complete graphs than repair
        // heuristics of ILP solvers
        /*
        if(i != 0) {
            visitor.repairSolution();
            repairSolution();
        }
        */
        visitor.optimize();
        ilp_.optimize();
        visitor.addCycleInequalities();
        std::size_t nc = addCycleInequalities();
        visitor.endIteration(i, nc);
        if(nc == 0) {
            break;
        }
    }

    visitor.repairSolution();
    repairSolution();
}

template<class GRAPH_VISITOR, class ILP>
inline typename Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::Ilp&
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::ilp() {
    return ilp_;
}

template<class GRAPH_VISITOR, class ILP>
inline double
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::label(
    const std::size_t edge
) const {
    assert(edge < graph_->numberOfEdges());
    size_t v0 = graph_->vertexOfEdge(edge, 0);
    size_t v1 = graph_->vertexOfEdge(edge, 1);
    return components_.areConnected(v0, v1) ? 0 : 1;
}

template<class GRAPH_VISITOR, class ILP>
inline std::size_t
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::addCycleInequalities() {
    std::size_t numberOfInequalitiesAdded = 0;
    std::size_t vi[] = {0, 0, 0};
    for(std::size_t j = 0; j < graph_->numberOfVertices(); ++j) {
        for(std::size_t k = j + 1; k < graph_->numberOfVertices(); ++k) {
            vi[0] = graph_->findEdge(j, k).second;
            for(std::size_t l = k + 1; l < graph_->numberOfVertices(); ++l) {
                vi[1] = graph_->findEdge(k, l).second;
                vi[2] = graph_->findEdge(j, l).second;
                const double lowerBound = 0.0;
                const double upperBound = std::numeric_limits<double>::infinity();
                if(ilp_.label(vi[0]) == 0) {
                    if(ilp_.label(vi[1]) == 0) {
                        if(ilp_.label(vi[2]) == 1) {
                            const double coefficients[] = {1.0, 1.0, -1.0};
                            ilp_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                            ++numberOfInequalitiesAdded;
                        }
                    }
                    else {
                        if(ilp_.label(vi[2]) == 0) {
                            const double coefficients[] = {1.0, -1.0, 1.0};
                            ilp_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                            ++numberOfInequalitiesAdded;
                        }
                    }
                }
                else {
                    if(ilp_.label(vi[1]) == 0 && ilp_.label(vi[2]) == 0) {
                        const double coefficients[] = {-1.0, 1.0, 1.0};
                        ilp_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                        ++numberOfInequalitiesAdded;
                    }
                }
            }
        }
    }
    return numberOfInequalitiesAdded;
}

template<class GRAPH_VISITOR, class ILP>
inline void
Multicut<CompleteGraph<GRAPH_VISITOR>, ILP>::repairSolution() {
    components_.build(*graph_, SubgraphWithoutCut(ilp_));
    std::vector<double> repairedSolution(graph_->numberOfEdges());
    for(std::size_t edge = 0; edge < graph_->numberOfEdges(); ++edge) {
        const std::size_t v0 = graph_->vertexOfEdge(edge, 0);
        const std::size_t v1 = graph_->vertexOfEdge(edge, 1);
        if(!components_.areConnected(v0, v1)) {
            repairedSolution[edge] = 1;
        }
    }
    ilp_.setStart(repairedSolution.begin());
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_HXX
