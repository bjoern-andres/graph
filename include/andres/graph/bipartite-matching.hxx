#pragma once
#ifndef ANDRES_GRAPH_MATCHING_HXX
#define ANDRES_GRAPH_MATCHING_HXX

#include <algorithm> // std::count
#include <cassert>
#include <cmath>  // std::abs
#include <limits> // std::numeric_limits::max() and ::epsilon()
#include <stdexcept>
#include <vector>

namespace andres {
namespace graph {

template <class GRAPH, class COST, class EDGELABELS>
inline void findMCBM(GRAPH const& graph, COST& costs, EDGELABELS& edge_labels);

/// Hungarian algorithm for the minimimum cost bipartite matching problem on a
/// Digraph.
///
/// Implementation based on:
/// Bourgeois and Lasalle, 1971, "An extension of the Munkres algorithm for the
/// assignment problem to rectangular matrices".
///
/// \tparam GRAPH A directed graph.
/// \tparam COST Edge cost container.
/// \tparam EDGELABEL Edge label container.
///
template <class GRAPH, class COST, class EDGELABELS>
class BipartiteMatching
{
public:
    static_assert(std::is_integral<typename EDGELABELS::value_type>::value,
                  "Integer-valued edge_labels required!");
    static_assert(!std::is_same<typename EDGELABELS::value_type, bool>::value,
                  "Edge_Labels cannot be boolean-valued!");

    struct IndexPair
    {
        IndexPair() {}
        IndexPair(std::size_t _first, std::size_t _second)
          : first(_first)
          , second(_second)
        {
        }
        typename std::size_t first, second;
    };

    BipartiteMatching(GRAPH const& graph, COST& costs, EDGELABELS& edge_labels);

    void run();
    std::vector<IndexPair> matches() const;

private:
    using value_type = typename COST::value_type;

    // members
    GRAPH const& graph_;
    COST& costs_;
    EDGELABELS& edge_labels_;

    std::vector<std::size_t> rows_;
    std::vector<std::size_t> cols_;
    std::vector<bool> row_cover_;
    std::vector<bool> col_cover_;
    std::vector<IndexPair> path_;
    std::size_t max_matches_;

    // labels
    constexpr static typename EDGELABELS::value_type NORMAL = 0;
    constexpr static typename EDGELABELS::value_type STARRED = 1;
    constexpr static typename EDGELABELS::value_type PRIMED = 2;

    enum class Step
    {
        ONE,
        TWO,
        THREE,
        FOUR,
        FIVE,
        SIX,
        SEVEN
    };

    Step step_{ Step::ONE };

    // basic steps of munkres algorithm.
    void step_one();
    void step_two();
    void step_three();
    void step_four();
    void step_five();
    void step_six();
    void step_seven();

    // auxiliary methods.
    IndexPair find_zero(std::size_t start_row = 0) const;
    bool has_star_in_row(std::size_t row) const;
    std::size_t find_star_in_row(std::size_t row) const;
    std::size_t find_star_in_col(std::size_t col) const;
    std::size_t find_prime_in_row(std::size_t row) const;
    value_type find_smallest() const;
    void clear_covers();
    void augment();
    void erase_primes();

    static constexpr std::size_t IDXMAX =
        std::numeric_limits<std::size_t>::max();
};

// implementation

/// \brief Find a minimimum cost bipartite matching by the hungarian algorithm.
///
/// Implementation based on:
/// Bourgeois and Lasalle, 1971, "An extension of the Munkres algorithm for the
/// assignment problem to rectangular matrices"
///
/// \param graph A bipartite, directed graph with all edges going from
/// one set of vertices to the other.
/// \param costs Edge costs. Will be modified in-place.
/// \param edge_labels Edge labels. Have to be initialized to 0 and will be
/// modified in-place.
///
template <class GRAPH, class COST, class EDGELABELS>
inline void
findMCBM(GRAPH const& graph, COST& costs, EDGELABELS& edge_labels)
{
    auto matcher =
        BipartiteMatching<GRAPH, COST, EDGELABELS>(graph, costs, edge_labels);
    matcher.run();
}

/// Construct an instance of the hungarian matching algorithm for biparite
/// Digraphs.
///
/// \param graph A bipartite, directed graph with all edges going from
/// one set of vertices to the other.
/// \param costs Edge costs. Will be modified in-place.
/// \param edge_labels Edge labels. Have to be initialized to 0 and will be
/// modified in-place.
///
template <class GRAPH, class COST, class EDGELABELS>
BipartiteMatching<GRAPH, COST, EDGELABELS>::BipartiteMatching(
    GRAPH const& graph, COST& costs, EDGELABELS& edge_labels)
  : graph_(graph)
  , costs_(costs)
  , edge_labels_(edge_labels)
{
    for (std::size_t v = 0; v < graph_.numberOfVertices(); ++v) {
        if (graph_.numberOfEdgesToVertex(v) > 0)
            cols_.emplace_back(v);
        else if (graph_.numberOfEdgesFromVertex(v) > 0)
            rows_.emplace_back(v);

        // debug only: check for bipartite-ness.
        assert(!(graph_.numberOfEdgesToVertex(v) > 0 &&
                 graph_.numberOfEdgesFromVertex(v) > 0));
    }

    row_cover_.resize(rows_.back() + 1, false);
    col_cover_.resize(cols_.back() + 1, false);

    max_matches_ = std::min(rows_.size(), cols_.size());
}

/// \brief Calculate matching.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::run()
{
    while (true) {
        switch (step_) {
            case Step::ONE:
                step_one();
                break;
            case Step::TWO:
                step_two();
                break;
            case Step::THREE:
                step_three();
                break;
            case Step::FOUR:
                step_four();
                break;
            case Step::FIVE:
                step_five();
                break;
            case Step::SIX:
                step_six();
                break;
            case Step::SEVEN:
                return;
            default:
                throw std::runtime_error("Failed to handle step code!");
        }
    }
}

/// \brief Generate pairs of matched vertices.
///
/// Requires BipartiteMatching::run() to be called first.
///
template <class GRAPH, class COST, class EDGELABELS>
std::vector<typename BipartiteMatching<GRAPH, COST, EDGELABELS>::IndexPair>
BipartiteMatching<GRAPH, COST, EDGELABELS>::matches() const
{
    std::vector<IndexPair> matches;
    for (auto const& row : rows_) {
        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it)
            if (edge_labels_[it->edge()] == STARRED)
                matches.emplace_back(row, it->vertex());
    }
    return matches;
}

/// subtract smallest value of each row/column.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::step_one()
{
    value_type minval;

    if (rows_.size() >= cols_.size()) {

        // subtract minval from each column.
        for (auto const& col : cols_) {

            minval = std::numeric_limits<value_type>::max();

            for (auto it = graph_.adjacenciesToVertexBegin(col);
                 it != graph_.adjacenciesToVertexEnd(col); ++it) {
                auto const& val = costs_[it->edge()];
                if (val < minval)
                    minval = val;
            }

            // subtract.
            for (auto it = graph_.adjacenciesToVertexBegin(col);
                 it != graph_.adjacenciesToVertexEnd(col); ++it) {
                costs_[it->edge()] -= minval;
            }
        }

    } else {

        // subtract minval from each row.
        for (auto const& row : rows_) {

            minval = std::numeric_limits<value_type>::max();

            for (auto it = graph_.adjacenciesFromVertexBegin(row);
                 it != graph_.adjacenciesFromVertexEnd(row); ++it) {
                auto const& val = costs_[it->edge()];
                if (val < minval)
                    minval = val;
            }

            for (auto it = graph_.adjacenciesFromVertexBegin(row);
                 it != graph_.adjacenciesFromVertexEnd(row); ++it) {
                costs_[it->edge()] -= minval;
            }
        }
    }
    step_ = Step::TWO;
}

/// Find zeros and star them.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::step_two()
{
    for (auto const& row : rows_) {
        if (row_cover_[row])
            continue;

        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it) {

            // star zero entries.
            if (!col_cover_[it->vertex()] && !row_cover_[row] &&
                costs_[it->edge()] == 0) {
                col_cover_[it->vertex()] = true;
                row_cover_[row] = true;
                edge_labels_[it->edge()] = STARRED;
            }
        }
    }

    // reset covers.
    clear_covers();

    step_ = Step::THREE;
}

/// Count covered columns and stop if k=min(n_rows, n_cols)
/// columns are covered.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::step_three()
{
    for (auto const& row : rows_)
        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it)
            if (edge_labels_[it->edge()] == STARRED)
                col_cover_[it->vertex()] = true;

    const std::size_t col_count =
        std::count(col_cover_.cbegin(), col_cover_.cend(), true);

    // if the matching is already perfect, then we are done.
    if (col_count >= max_matches_)
        step_ = Step::SEVEN;
    else
        step_ = Step::FOUR;
}

/// Find non-covered zeros and prime them.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::step_four()
{
    IndexPair idx;
    std::size_t start_row{ 0 };

    while (true) {
        idx = find_zero(start_row);

        if (idx.first == IDXMAX || idx.second == IDXMAX) {
            // No zero found.
            step_ = Step::SIX;
            return;
        }

        // it is more likely to find another zero
        // _after_ the last location, so we pass it
        // as a hint to find_zero.
        start_row = idx.first + 1;

        // prime at given index.
        edge_labels_[graph_.findEdge(idx.first, idx.second).second] = PRIMED;

        if (has_star_in_row(idx.first)) {
            const auto col = find_star_in_row(idx.first);
            row_cover_[idx.first] = true;
            col_cover_[col] = false;
        } else {
            // starting point for augmenting path.
            path_.clear();
            path_.emplace_back(idx.first, idx.second);
            step_ = Step::FIVE;
            return;
        }
    }
}

/// construct augmenting path of starred and primed zeros.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::step_five()
{
    while (true) {
        auto const row = find_star_in_col(path_.back().second);
        if (row == IDXMAX) {
            break;
        }

        path_.emplace_back(row, path_.back().second);
        path_.emplace_back(row, find_prime_in_row(row));
    }

    augment();
    clear_covers();
    erase_primes();
    step_ = Step::THREE;
}

/// update cost based on non-satisfied constraints.
///
template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::step_six()
{
    auto const& minval = find_smallest();
    for (auto const& row : rows_) {
        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it) {
            if (row_cover_[row])
                costs_[it->edge()] += minval;
            if (!col_cover_[it->vertex()])
                costs_[it->edge()] -= minval;
        }
    }
    step_ = Step::FOUR;
}

template <class GRAPH, class COST, class EDGELABELS>
typename BipartiteMatching<GRAPH, COST, EDGELABELS>::value_type
BipartiteMatching<GRAPH, COST, EDGELABELS>::find_smallest() const
{
    auto minval = std::numeric_limits<value_type>::max();
    for (auto const& row : rows_) {
        // consider only non-covered rows.
        if (row_cover_[row])
            continue;
        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it)
            // consider only non-covered columns.
            if (!col_cover_[it->vertex()] && costs_[it->edge()] < minval)
                minval = costs_[it->edge()];
    }
    return minval;
}

template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::augment()
{
    for (auto const& idx : path_) {
        auto const& edge = graph_.findEdge(idx.first, idx.second).second;
        if (edge_labels_[edge] == STARRED)
            edge_labels_[edge] = NORMAL;
        else
            edge_labels_[edge] = STARRED;
    }
}

template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::clear_covers()
{
    for (std::size_t row = 0; row < row_cover_.size(); ++row) {
        row_cover_[row] = false;
    }
    for (std::size_t col = 0; col < col_cover_.size(); ++col) {
        col_cover_[col] = false;
    }
}

template <class GRAPH, class COST, class EDGELABELS>
void
BipartiteMatching<GRAPH, COST, EDGELABELS>::erase_primes()
{
    for (auto const& row : rows_)
        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it)
            if (edge_labels_[it->edge()] == PRIMED)
                edge_labels_[it->edge()] = NORMAL;
}

template <class GRAPH, class COST, class EDGELABELS>
typename BipartiteMatching<GRAPH, COST, EDGELABELS>::IndexPair
BipartiteMatching<GRAPH, COST, EDGELABELS>::find_zero(
    std::size_t start_row) const
{
    auto const& n_rows = rows_.size();

    std::size_t row_idx = start_row;
    for (std::size_t idx = 0; idx < n_rows; ++idx, ++row_idx) {

        // start from start_row, and wrap around at n_rows.
        if (row_idx >= n_rows) {
            row_idx = 0;
        }

        auto const& row = rows_[row_idx];

        if (row_cover_[row])
            continue;

        for (auto it = graph_.adjacenciesFromVertexBegin(row);
             it != graph_.adjacenciesFromVertexEnd(row); ++it) {
            if (!col_cover_[it->vertex()] &&
                std::abs(costs_[it->edge()]) <=
                    std::numeric_limits<value_type>::epsilon())
                return { row, it->vertex() };
        }
    }

    return { IDXMAX, IDXMAX };
}

template <class GRAPH, class COST, class EDGELABELS>
bool
BipartiteMatching<GRAPH, COST, EDGELABELS>::has_star_in_row(
    std::size_t const row) const
{
    for (auto it = graph_.adjacenciesFromVertexBegin(row);
         it != graph_.adjacenciesFromVertexEnd(row); ++it)
        if (edge_labels_[it->edge()] == STARRED)
            return true;
    return false;
}

template <class GRAPH, class COST, class EDGELABELS>
std::size_t
BipartiteMatching<GRAPH, COST, EDGELABELS>::find_star_in_row(
    std::size_t const row) const
{
    for (auto it = graph_.adjacenciesFromVertexBegin(row);
         it != graph_.adjacenciesFromVertexEnd(row); ++it)
        if (edge_labels_[it->edge()] == STARRED)
            return it->vertex();

    return IDXMAX;
}

template <class GRAPH, class COST, class EDGELABELS>
std::size_t
BipartiteMatching<GRAPH, COST, EDGELABELS>::find_star_in_col(
    std::size_t const col) const
{
    for (auto it = graph_.adjacenciesToVertexBegin(col);
         it != graph_.adjacenciesToVertexEnd(col); ++it)
        if (edge_labels_[it->edge()] == STARRED)
            return it->vertex();

    return IDXMAX;
}

template <class GRAPH, class COST, class EDGELABELS>
std::size_t
BipartiteMatching<GRAPH, COST, EDGELABELS>::find_prime_in_row(
    std::size_t const row) const
{
    for (auto it = graph_.adjacenciesFromVertexBegin(row);
         it != graph_.adjacenciesFromVertexEnd(row); ++it)
        if (edge_labels_[it->edge()] == PRIMED)
            return it->vertex();

    throw std::runtime_error("Could not find prime in given row!");
}

} // end namespace graph
} // end namespace andres

#endif
