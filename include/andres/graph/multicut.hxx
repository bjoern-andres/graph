#pragma once
#ifndef ANDRES_MULTICUT_HXX
#define ANDRES_MULTICUT_HXX

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <deque>
#include <algorithm> // std::copy

#include "andres/partition.hxx" 
#include "andres/graph/paths.hxx"
#include "andres/graph/components.hxx"
#include "andres/graph/shortest-paths.hxx"

namespace andres {
namespace graph {

/// Multicut solver.
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
template<class GRAPH, class ILP, class COMPONENTS = ComponentsBySearch<GRAPH> >
class Multicut {
public:
    typedef GRAPH Graph;
    typedef ILP Ilp;
    typedef COMPONENTS Components;
    
    /// Visitors can be used to follow the progress of the optimization.
    struct IdleVisitor {
        void repairSolution() const {}
        void optimize() const {}
        void addCycleInequalities() const {}
        void endIteration(const std::size_t, const std::size_t) const {}
    };

    Multicut();
    void setup(const Graph&, const std::vector<double>&);
    void solve(const std::size_t);
    template<class VISITOR>
        void solve(const std::size_t, VISITOR&);
    
    Ilp& ilp();
    double label(const std::size_t) const;

private:
    struct SubgraphWithoutCut { // a subgraph mask
        SubgraphWithoutCut(const Multicut<GRAPH, ILP, COMPONENTS>& multicut)
            : multicut_(multicut) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            { return multicut_.label(e) == 0; }

        const Multicut<GRAPH, ILP, COMPONENTS>& multicut_;
    };

    std::size_t addCycleInequalities();
    void repairSolution();

    Components components_;
    const Graph* graph_;
    Ilp ilp_;
};

template<class GRAPH, class ILP, class COMPONENTS>
inline
Multicut<GRAPH, ILP, COMPONENTS>::Multicut() 
:   components_(),
    graph_(NULL),
    ilp_()
{}

template<class GRAPH, class ILP, class COMPONENTS>
inline double
Multicut<GRAPH, ILP, COMPONENTS>::label(
    const std::size_t edge
) const {
    assert(edge < graph_->numberOfEdges());
    return ilp_.label(edge);
}

template<class GRAPH, class ILP, class COMPONENTS>
inline typename Multicut<GRAPH, ILP, COMPONENTS>::Ilp&
Multicut<GRAPH, ILP, COMPONENTS>::ilp() {
    return ilp_;
}

template<class GRAPH, class ILP, class COMPONENTS>
void 
Multicut<GRAPH, ILP, COMPONENTS>::setup(
    const Graph& graph,
    const std::vector<double>& weights    
) {
    if(graph.numberOfEdges() != weights.size()) {
        throw std::runtime_error("the number of edges in the graph does not equal the number of weights.");
    }
    graph_ = &graph;
    ilp_.initModel(graph_->numberOfEdges(), weights.data());
}

template<class GRAPH, class ILP, class COMPONENTS>
inline void 
Multicut<GRAPH, ILP, COMPONENTS>::solve(
    const std::size_t maxIterations
) {
    solve(maxIterations, IdleVisitor());
}

template<class GRAPH, class ILP, class COMPONENTS>
template<class VISITOR>
void 
Multicut<GRAPH, ILP, COMPONENTS>::solve(
    const std::size_t maxIterations,
    VISITOR& visitor
) {
    for(std::size_t i = 0; maxIterations == 0 || i < maxIterations; ++i) {
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
}

template<class GRAPH, class ILP, class COMPONENTS>
std::size_t
Multicut<GRAPH, ILP, COMPONENTS>::addCycleInequalities() {
    // label connected components
    components_.build(*graph_, SubgraphWithoutCut(*this));

    // search for violated non-chordal cycles and add corresp. inequalities
    std::deque<std::size_t> path;
    std::vector<ptrdiff_t> buffer;
    std::vector<double> variables;
    std::vector<double> coefficients;
    std::size_t counter = 0;
    #pragma omp parallel for private(path, buffer, variables, coefficients), schedule(guided)
    for(ptrdiff_t edge = 0; edge < graph_->numberOfEdges(); ++edge) {
        if(label(edge) == 1) {
            const std::size_t v0 = graph_->vertexOfEdge(edge, 0);
            const std::size_t v1 = graph_->vertexOfEdge(edge, 1);
            if(components_.areConnected(v0, v1)) { 
                // search for shortest path
                bool found = spsp(*graph_, SubgraphWithoutCut(*this), v0, v1, path, buffer); 
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
template<class GRAPH, class ILP, class COMPONENTS>
void
Multicut<GRAPH, ILP, COMPONENTS>::repairSolution() {
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

#endif // #ifndef ANDRES_MULTICUT_HXX
