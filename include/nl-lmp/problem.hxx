#pragma once
#ifndef NL_LMP_PROBLEM_HXX
#define NL_LMP_PROBLEM_HXX

#include <andres/marray.hxx>
#include <andres/graph/complete-graph.hxx>



namespace nl_lmp
{

template<typename GRAPH>
class Problem
{
public:
    typedef GRAPH GraphType;
    typedef GRAPH LiftedGraphType;


    Problem(size_t numberOfVertices, size_t numberOfClasses)
    {
        unaryCosts_.resize({ numberOfVertices, numberOfClasses });

        original_graph_.assign(numberOfVertices);
        lifted_graph_.assign(numberOfVertices);
    }

    double getUnaryCost(size_t v, size_t c) const
    {
        return unaryCosts_(v, c);
    }

    double getPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseCutCosts_[liftedGraph().findEdge(v0, v1).second](c0, c1);
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    double getPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1, size_t edge_index) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseCutCosts_[edge_index](c0, c1);
    }

    double getPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseJoinCosts_[liftedGraph().findEdge(v0, v1).second](c0, c1);
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    double getPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1, size_t edge_index) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseJoinCosts_[edge_index](c0, c1);
    }

    GraphType const& originalGraph() const
    {
        return original_graph_;
    }

    LiftedGraphType const& liftedGraph() const
    {
        return lifted_graph_;
    }

    size_t numberOfClasses() const
    {
        return unaryCosts_.shape(1);
    }

    size_t numberOfVertices() const
    {
        return unaryCosts_.shape(0);
    }

    void setUnaryCost(size_t v, size_t c, double value)
    {
        unaryCosts_(v, c) = value;
    }

    void setPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1, double value, bool add_edge_into_original_graph = true)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        if (add_edge_into_original_graph)
            original_graph_.insertEdge(v0, v1);

        auto edge_index = lifted_graph_.insertEdge(v0, v1);
        if (edge_index == pairwiseCutCosts_.size())
        {
            pairwiseCutCosts_.push_back(andres::Marray<double>({ numberOfClasses(), numberOfClasses() }));
            pairwiseJoinCosts_.push_back(andres::Marray<double>({ numberOfClasses(), numberOfClasses() }));
        }

        pairwiseCutCosts_[edge_index](c0, c1) = value;
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    void setPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1, double value, size_t edge_index)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        pairwiseCutCosts_[edge_index](c0, c1) = value;
    }

    void setPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1, double value, bool add_edge_into_original_graph = true)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        if (add_edge_into_original_graph)
            original_graph_.insertEdge(v0, v1);

        auto edge_index = lifted_graph_.insertEdge(v0, v1);
        if (edge_index == pairwiseCutCosts_.size())
        {
            pairwiseCutCosts_.push_back(andres::Marray<double>({ numberOfClasses(), numberOfClasses() }));
            pairwiseJoinCosts_.push_back(andres::Marray<double>({ numberOfClasses(), numberOfClasses() }));
        }

        pairwiseJoinCosts_[edge_index](c0, c1) = value;
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    void setPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1, double value, size_t edge_index)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        pairwiseJoinCosts_[edge_index](c0, c1) = value;
    }

private:
    GraphType original_graph_;
    LiftedGraphType lifted_graph_;

    std::vector<andres::Marray<double>> pairwiseCutCosts_;
    std::vector<andres::Marray<double>> pairwiseJoinCosts_;

    andres::Marray<double> unaryCosts_;
};

template<typename GRAPHVISITOR>
class Problem<andres::graph::CompleteGraph<GRAPHVISITOR>>
{
public:
    typedef andres::graph::CompleteGraph<GRAPHVISITOR> GraphType;


    Problem(size_t numberOfVertices, size_t numberOfClasses)
    {
        graph_.assign(numberOfVertices);

        unaryCosts_.resize({ graph_.numberOfVertices(), numberOfClasses });

        pairwiseCutCosts_.resize(graph_.numberOfEdges());
        pairwiseJoinCosts_.resize(graph_.numberOfEdges());

        for (size_t e = 0; e < graph_.numberOfEdges(); ++e)
        {
            pairwiseCutCosts_[e].resize({ numberOfClasses, numberOfClasses });
            pairwiseJoinCosts_[e].resize({ numberOfClasses, numberOfClasses });
        }
    }

    double getUnaryCost(size_t v, size_t c) const
    {
        return unaryCosts_(v, c);
    }

    double getPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseCutCosts_[graph_.findEdge(v0, v1).second](c0, c1);
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    double getPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1, size_t edge_index) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseCutCosts_[edge_index](c0, c1);
    }

    double getPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseJoinCosts_[graph_.findEdge(v0, v1).second](c0, c1);
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    double getPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1, size_t edge_index) const
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return pairwiseJoinCosts_[edge_index](c0, c1);
    }

    GraphType const& originalGraph() const
    {
        return graph_;
    }

    GraphType const& liftedGraph() const
    {
        return graph_;
    }

    size_t numberOfClasses() const
    {
        return unaryCosts_.shape(1);
    }

    size_t numberOfEdges() const
    {
        return pairwiseCutCosts_.size();
    }

    size_t numberOfVertices() const
    {
        return unaryCosts_.shape(0);
    }

    void setUnaryCost(size_t v, size_t c, double value)
    {
        unaryCosts_(v, c) = value;
    }

    void setPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1, double value)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        pairwiseCutCosts_[graph_.findEdge(v0, v1).second](c0, c1) = value;
    }

    // by using this method you agree to all consequences due to possible wrong indexing...
    void setPairwiseCutCost(size_t v0, size_t v1, size_t c0, size_t c1, double value, size_t edge_index)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        pairwiseCutCosts_[edge_index](c0, c1) = value;
    }

    // for complete graph these methods do nothing
    // they are necessary only for unified interface, so far...
    void setPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1, double value)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        pairwiseJoinCosts_[graph_.findEdge(v0, v1).second](c0, c1) = value;
    }

    void setPairwiseJoinCost(size_t v0, size_t v1, size_t c0, size_t c1, double value, size_t edge_index)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        pairwiseJoinCosts_[edge_index](c0, c1) = value;
    }

private:
    GraphType graph_;

    std::vector<andres::Marray<double>> pairwiseCutCosts_;
    std::vector<andres::Marray<double>> pairwiseJoinCosts_;
    andres::Marray<double> unaryCosts_;
};

}

#endif
