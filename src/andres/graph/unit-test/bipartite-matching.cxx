#include <stdexcept>
#include <utility> // std::pair
#include <vector>

#include "andres/graph/bipartite-matching.hxx"
#include "andres/graph/digraph.hxx"

inline void
test(const bool& pred)
{
    if (!pred)
        throw std::runtime_error("Test failed.");
}

template <class GRAPH, class COST, class MATCHES>
void
test_case(GRAPH& graph, COST& cost, MATCHES const& expected_matches)
{
    using labels_t = std::vector<unsigned int>;
    using pair_t = std::pair<std::size_t, std::size_t>;

    // copy costs to calculate total sum of matchings.
    const auto original_costs = cost;

    // find best matching.
    labels_t edge_labels(cost.size(), 0);
    andres::graph::findMCBM(graph, cost, edge_labels);

    // and construct matches from mask.
    const auto matches = [](GRAPH const& graph, labels_t const& edge_labels) {
        std::vector<pair_t> matches;
        for (std::size_t edge = 0; edge < edge_labels.size(); ++edge) {
            if (edge_labels[edge] == 1) {
                const auto v0 = graph.vertexOfEdge(edge, 0);
                const auto v1 = graph.vertexOfEdge(edge, 1);
                matches.emplace_back(std::move(v0), std::move(v1));
            }
        }
        return matches;
    }(graph, edge_labels);

    // compare matching with expected outcome.
    typename COST::value_type sum = 0;
    for (std::size_t idx = 0; idx < cost.size(); ++idx)
        if (edge_labels[idx] == 1)
            sum += original_costs[idx];

    typename COST::value_type expected_sum = 0;
    for (auto const& match : expected_matches)
        expected_sum +=
            original_costs[graph.findEdge(match.first, match.second).second];

    test(sum == expected_sum);

    for (auto const& match : matches) {
        auto found = std::find_if(
            expected_matches.cbegin(), expected_matches.cend(),
            [=](pair_t const& pair) {
                if (pair.second != match.second || pair.first != match.first)
                    return false;
                else
                    return true;
            });
        test(found != expected_matches.cend());
    }

    test(matches.size() == expected_matches.size());
}

int
main()
{
    {
        using cost_t = std::vector<double>;
        using graph_t = andres::graph::Digraph<>;
        using pair_t = std::pair<std::size_t, std::size_t>;

        // construct graph and costs.
        constexpr std::size_t n_rows = 4;
        constexpr std::size_t n_cols = 6;

        graph_t graph(n_rows + n_cols);
        cost_t costs;
        costs.reserve(n_rows * n_cols);

        for (std::size_t row = 0; row < n_rows; ++row) {
            for (std::size_t col = 0; col < n_cols; ++col) {
                graph.insertEdge(row, col + n_rows);
                costs.emplace_back((row + 1) * (col + 1));
            }
        }

        // and the expected matches.
        const std::vector<pair_t> expected_matches{ { 0, 3 + n_rows },
                                                    { 1, 2 + n_rows },
                                                    { 2, 1 + n_rows },
                                                    { 3, 0 + n_rows } };

        test_case(graph, costs, expected_matches);
    }

    {
        using cost_t = std::vector<double>;
        using graph_t = andres::graph::Digraph<>;
        using pair_t = std::pair<std::size_t, std::size_t>;

        // construct graph and costs.
        graph_t graph(9);
        cost_t costs;
        costs.reserve(7);

        auto add_edge = [&](std::size_t v, std::size_t w,
                            typename cost_t::value_type cost) {
            graph.insertEdge(v, w);
            costs.emplace_back(std::move(cost));
        };

        // no edges for v=0.
        add_edge(1, 4, 3.0); // x
        add_edge(1, 5, 5.0);
        add_edge(2, 5, 5.0);
        add_edge(2, 6, -1.0); // x
        add_edge(3, 6, -2.0);
        add_edge(3, 7, .0);
        add_edge(3, 8, -1.0); // x

        // and the expected matches.
        const std::vector<pair_t> expected_matches{ { 1, 4 },
                                                    { 2, 6 },
                                                    { 3, 8 } };

        test_case(graph, costs, expected_matches);
    }

    {
        using cost_t = std::vector<double>;
        using graph_t = andres::graph::Digraph<>;
        using pair_t = std::pair<std::size_t, std::size_t>;

        // construct graph and costs.
        graph_t graph(9);
        cost_t costs;
        costs.reserve(7);

        auto add_edge = [&](std::size_t v, std::size_t w,
                            typename cost_t::value_type cost) {
            graph.insertEdge(v, w);
            costs.emplace_back(std::move(cost));
        };

        // no edges for v=0.
        add_edge(4, 1, 3.0); // x
        add_edge(5, 1, 5.0);
        add_edge(5, 2, 5.0);
        add_edge(6, 2, -1.0); // x
        add_edge(6, 3, -2.0);
        add_edge(7, 3, .0);
        add_edge(8, 3, -1.0); // x

        // and the expected matches.
        const std::vector<pair_t> expected_matches{ { 4, 1 },
                                                    { 6, 2 },
                                                    { 8, 3 } };

        test_case(graph, costs, expected_matches);
    }

    {
        using cost_t = std::vector<double>;
        using graph_t = andres::graph::Digraph<>;
        using pair_t = std::pair<std::size_t, std::size_t>;

        // construct graph and costs.

        graph_t graph(30);
        cost_t costs;

        auto add_edge = [&](std::size_t v, std::size_t w,
                            typename cost_t::value_type cost) {
            graph.insertEdge(v, w);
            costs.emplace_back(std::move(cost));
        };

        add_edge(0, 15, 0.48);
        add_edge(0, 16, 1.2);
        add_edge(0, 17, 0.38);
        add_edge(0, 18, -0.07);
        add_edge(0, 19, -1.07);
        add_edge(0, 20, 0.68);
        add_edge(0, 21, -0.96);
        add_edge(0, 22, 1.09);
        add_edge(0, 23, 1.69);
        add_edge(0, 24, -0.68);
        add_edge(0, 25, 1.37);
        add_edge(0, 26, -0.85);
        add_edge(0, 27, -1.14);
        add_edge(0, 28, 1.21);
        add_edge(0, 29, -1.48);
        add_edge(1, 15, -0.75);
        add_edge(1, 16, 0.32);
        add_edge(1, 17, -0.54);
        add_edge(1, 18, 0.11);
        add_edge(1, 19, -1.41);
        add_edge(1, 20, 0.69);
        add_edge(1, 21, -1.53);
        add_edge(1, 22, -0.42);
        add_edge(1, 23, 0.13);
        add_edge(1, 24, -0.44);
        add_edge(1, 25, -1.29);
        add_edge(1, 26, 1.9);
        add_edge(1, 27, 0.92);
        add_edge(1, 28, -1.11);
        add_edge(1, 29, -0.29);
        add_edge(2, 15, 2.29);
        add_edge(2, 16, 0.25);
        add_edge(2, 17, 0.36);
        add_edge(2, 18, 2.83);
        add_edge(2, 19, -0.5);
        add_edge(2, 20, -0.42);
        add_edge(2, 21, -0.1);
        add_edge(2, 22, -1.23);
        add_edge(2, 23, 0.32);
        add_edge(2, 24, -0.22);
        add_edge(2, 25, -0.15);
        add_edge(2, 26, 0.35);
        add_edge(2, 27, 2.38);
        add_edge(2, 28, -0.93);
        add_edge(2, 29, 1.28);
        add_edge(3, 15, -1.06);
        add_edge(3, 16, 0.67);
        add_edge(3, 17, 1.1);
        add_edge(3, 18, 0.77);
        add_edge(3, 19, -1.05);
        add_edge(3, 20, 0.36);
        add_edge(3, 21, -1.09);
        add_edge(3, 22, -1.77);
        add_edge(3, 23, -0.35);
        add_edge(3, 24, -1.13);
        add_edge(3, 25, 0.8);
        add_edge(3, 26, -1.96);
        add_edge(3, 27, -0.5);
        add_edge(3, 28, 0.9);
        add_edge(3, 29, -0.13);
        add_edge(4, 15, 0.73);
        add_edge(4, 16, -0.28);
        add_edge(4, 17, 0.86);
        add_edge(4, 18, 1.36);
        add_edge(4, 19, 0.84);
        add_edge(4, 20, -0.01);
        add_edge(4, 21, -0.6);
        add_edge(4, 22, 1.09);
        add_edge(4, 23, -0.12);
        add_edge(4, 24, -0.18);
        add_edge(4, 25, 1.12);
        add_edge(4, 26, -0.99);
        add_edge(4, 27, 0.16);
        add_edge(4, 28, 1.44);
        add_edge(4, 29, -0.82);
        add_edge(5, 15, -0.45);
        add_edge(5, 16, -1.33);
        add_edge(5, 17, 0.32);
        add_edge(5, 18, 0.41);
        add_edge(5, 19, 0.82);
        add_edge(5, 20, 0.09);
        add_edge(5, 21, 2.18);
        add_edge(5, 22, -0.05);
        add_edge(5, 23, -0.17);
        add_edge(5, 24, 0.95);
        add_edge(5, 25, -2.26);
        add_edge(5, 26, 0.25);
        add_edge(5, 27, -2.02);
        add_edge(5, 28, 0.04);
        add_edge(5, 29, -0.39);
        add_edge(6, 15, -1.16);
        add_edge(6, 16, 0.59);
        add_edge(6, 17, -1.05);
        add_edge(6, 18, 2.01);
        add_edge(6, 19, -0.51);
        add_edge(6, 20, 0.3);
        add_edge(6, 21, -0.32);
        add_edge(6, 22, 1.03);
        add_edge(6, 23, 0.19);
        add_edge(6, 24, 0.41);
        add_edge(6, 25, -0.37);
        add_edge(6, 26, -1.01);
        add_edge(6, 27, -0.43);
        add_edge(6, 28, 0.2);
        add_edge(6, 29, 0.24);
        add_edge(7, 15, 0.48);
        add_edge(7, 16, -0.36);
        add_edge(7, 17, -1.59);
        add_edge(7, 18, 1.38);
        add_edge(7, 19, 0.73);
        add_edge(7, 20, -1.09);
        add_edge(7, 21, -1.4);
        add_edge(7, 22, 0.59);
        add_edge(7, 23, 0.76);
        add_edge(7, 24, 0.5);
        add_edge(7, 25, 1.08);
        add_edge(7, 26, -0.18);
        add_edge(7, 27, 0.72);
        add_edge(7, 28, 0.21);
        add_edge(7, 29, -0.88);
        add_edge(8, 15, 1.87);
        add_edge(8, 16, 0.59);
        add_edge(8, 17, 0.26);
        add_edge(8, 18, -0.44);
        add_edge(8, 19, 0.42);
        add_edge(8, 20, 0.01);
        add_edge(8, 21, 0.58);
        add_edge(8, 22, -0.53);
        add_edge(8, 23, 1.34);
        add_edge(8, 24, 0.96);
        add_edge(8, 25, 0.44);
        add_edge(8, 26, 2.31);
        add_edge(8, 27, 0.18);
        add_edge(8, 28, -0.89);
        add_edge(8, 29, -0.04);
        add_edge(9, 15, 1.06);
        add_edge(9, 16, 1.29);
        add_edge(9, 17, 0.03);
        add_edge(9, 18, -1.2);
        add_edge(9, 19, 0.65);
        add_edge(9, 20, -1.79);
        add_edge(9, 21, -0.36);
        add_edge(9, 22, -0.5);
        add_edge(9, 23, -0.11);
        add_edge(9, 24, -1.5);
        add_edge(9, 25, 1.27);
        add_edge(9, 26, -0.42);
        add_edge(9, 27, -0.52);
        add_edge(9, 28, -0.24);
        add_edge(9, 29, 0.15);
        add_edge(10, 15, 0.15);
        add_edge(10, 16, -1.21);
        add_edge(10, 17, -0.22);
        add_edge(10, 18, 0.82);
        add_edge(10, 19, -0.57);
        add_edge(10, 20, -0.33);
        add_edge(10, 21, -0.22);
        add_edge(10, 22, 0.75);
        add_edge(10, 23, 0.59);
        add_edge(10, 24, -1.95);
        add_edge(10, 25, -0.82);
        add_edge(10, 26, 0.32);
        add_edge(10, 27, 0.87);
        add_edge(10, 28, 0.33);
        add_edge(10, 29, 0.28);
        add_edge(11, 15, 1.31);
        add_edge(11, 16, -0.93);
        add_edge(11, 17, 0.17);
        add_edge(11, 18, -0.48);
        add_edge(11, 19, 0.33);
        add_edge(11, 20, 1.58);
        add_edge(11, 21, -0.75);
        add_edge(11, 22, -0.75);
        add_edge(11, 23, 1.56);
        add_edge(11, 24, -0.31);
        add_edge(11, 25, -1.44);
        add_edge(11, 26, 0.45);
        add_edge(11, 27, 0.55);
        add_edge(11, 28, -0.22);
        add_edge(11, 29, -0.86);
        add_edge(12, 15, 0.44);
        add_edge(12, 16, -1.11);
        add_edge(12, 17, 0.55);
        add_edge(12, 18, -1.54);
        add_edge(12, 19, -1.15);
        add_edge(12, 20, 0.78);
        add_edge(12, 21, -0.19);
        add_edge(12, 22, -1.33);
        add_edge(12, 23, 0.55);
        add_edge(12, 24, 0.78);
        add_edge(12, 25, -0.73);
        add_edge(12, 26, 0.51);
        add_edge(12, 27, -2.63);
        add_edge(12, 28, 0.73);
        add_edge(12, 29, -0.31);
        add_edge(13, 15, 0.47);
        add_edge(13, 16, 0.68);
        add_edge(13, 17, -0.98);
        add_edge(13, 18, 0.31);
        add_edge(13, 19, -0.65);
        add_edge(13, 20, 0.93);
        add_edge(13, 21, 0.16);
        add_edge(13, 22, -0.2);
        add_edge(13, 23, -1.3);
        add_edge(13, 24, -0.07);
        add_edge(13, 25, 1.91);
        add_edge(13, 26, 1.43);
        add_edge(13, 27, 0.32);
        add_edge(13, 28, -0.06);
        add_edge(13, 29, 0.01);
        add_edge(14, 15, -0.46);
        add_edge(14, 16, -1.34);
        add_edge(14, 17, -0.94);
        add_edge(14, 18, -0.54);
        add_edge(14, 19, -0.8);
        add_edge(14, 20, -1.03);
        add_edge(14, 21, -0.82);
        add_edge(14, 22, 0.27);
        add_edge(14, 23, 1.76);
        add_edge(14, 24, -2.11);
        add_edge(14, 25, -0.67);
        add_edge(14, 26, -0.45);
        add_edge(14, 27, 0.78);
        add_edge(14, 28, -1.43);
        add_edge(14, 29, 0.48);

        // and the expected matches.
        // Note: these reference matches were computed with pythons munkres
        // package.
        const std::vector<pair_t> expected_matches{
            { 0, 29 },  { 1, 19 },  { 2, 22 },  { 3, 26 },  { 4, 21 },
            { 5, 25 },  { 6, 15 },  { 7, 17 },  { 8, 18 },  { 9, 20 },
            { 10, 24 }, { 11, 16 }, { 12, 27 }, { 13, 23 }, { 14, 28 },
        };
        test_case(graph, costs, expected_matches);
    }

    return 0;
}
