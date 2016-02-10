#pragma once
#ifndef ANDRES_GRAPH_TWOCUT_LIFTED_KERNIGHAN_LIN_HXX
#define ANDRES_GRAPH_TWOCUT_LIFTED_KERNIGHAN_LIN_HXX

#include <limits>
#include <vector>

namespace andres {
namespace graph {
namespace twocut_lifted {

struct TwoCutSettings {
    std::size_t numberOfIterations { std::numeric_limits<std::size_t>::max() };
    double epsilon { 1e-9 };
};

template<typename ORIGINAL_GRAPH>
struct TwoCutBuffers {
    TwoCutBuffers(const ORIGINAL_GRAPH& graph)
        :   differences(graph.numberOfVertices()),
            is_moved(graph.numberOfVertices()),
            referenced_by(graph.numberOfVertices()),
            vertex_labels(graph.numberOfVertices())
        {}
    std::vector<double> differences;
    std::vector<char> is_moved;
    std::size_t max_not_used_label;
    std::vector<std::size_t> referenced_by;
    std::vector<std::size_t> vertex_labels;
};

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename SET, typename ECA>
inline double
kernighanLin(
    const ORIGINAL_GRAPH& original_graph,
    const LIFTED_GRAPH& lifted_graph,
    const ECA& edge_costs,
    SET& A,
    SET& B,
    TwoCutBuffers<ORIGINAL_GRAPH>& buffer,
    const TwoCutSettings settings
)
{
    struct Move
    {
        int v { -1 };
        double difference { std::numeric_limits<double>::lowest() };
        std::size_t new_label;
    };

    auto gain_from_merging = .0;

    auto compute_differences = [&](const SET& A, std::size_t label_A, std::size_t label_B)
    {
        for (long int i = 0; i < A.size(); ++i)
        {
            double diffExt = .0;
            double diffInt = .0;
            std::size_t ref_cnt = 0;

            for (auto it = lifted_graph.adjacenciesFromVertexBegin(A[i]); it != lifted_graph.adjacenciesFromVertexEnd(A[i]); ++it)
            {
                const auto lbl = buffer.vertex_labels[it->vertex()];

                if (lbl == label_A)
                    diffInt += edge_costs[it->edge()];
                else if (lbl == label_B)
                    diffExt += edge_costs[it->edge()];
            }

            for (auto it = original_graph.adjacenciesFromVertexBegin(A[i]); it != original_graph.adjacenciesFromVertexEnd(A[i]); ++it)
                if (buffer.vertex_labels[it->vertex()] == label_B)
                    ++ref_cnt;

            buffer.differences[A[i]] = diffExt - diffInt;
            buffer.referenced_by[A[i]] = ref_cnt;
            buffer.is_moved[A[i]] = 0;

            gain_from_merging += diffExt;
        }
    };


    if (A.empty())
        return .0;
    
    auto label_A = buffer.vertex_labels[A[0]];
    auto label_B = (!B.empty()) ? buffer.vertex_labels[B[0]] : buffer.max_not_used_label;

    compute_differences(A, label_A, label_B);
    compute_differences(B, label_B, label_A);

    gain_from_merging /= 2.0;

    std::vector<std::size_t> border;

    for (auto a : A)
        if (buffer.referenced_by[a] > 0)
            border.push_back(a);

    for (auto b : B)
        if (buffer.referenced_by[b] > 0)
            border.push_back(b);

    std::vector<Move> moves;
    double cumulative_diff = .0;
    std::pair<double, std::size_t> max_move { std::numeric_limits<double>::lowest(), 0 };

    for (std::size_t k = 0; k < settings.numberOfIterations; ++k)
    {
        Move m;

        if (B.empty() && k == 0)
        {
            for (auto a : A)
                if (buffer.differences[a] > m.difference)
                {
                    m.v = a;
                    m.difference = buffer.differences[a];
                }
        }
        else
        {
            std::size_t size = border.size();
            
            for (std::size_t i = 0; i < size; )
                if (buffer.referenced_by[border[i]] == 0)
                    std::swap(border[i], border[--size]);
                else
                {
                    if (buffer.differences[border[i]] > m.difference)
                    {
                        m.v = border[i];
                        m.difference = buffer.differences[m.v];
                    }
                    
                    ++i;
                }

            border.erase(border.begin() + size, border.end());
        }

        if (m.v == -1)
            break;

        const auto old_label = buffer.vertex_labels[m.v];

        if (old_label == label_A)
            m.new_label = label_B;
        else
            m.new_label = label_A;

        // update differences and references
        for (auto it = lifted_graph.adjacenciesFromVertexBegin(m.v); it != lifted_graph.adjacenciesFromVertexEnd(m.v); ++it)
        {
            if (buffer.is_moved[it->vertex()])
                continue;

            const auto lbl = buffer.vertex_labels[it->vertex()];
            
            // edge to an element of the new set
            if (lbl == m.new_label)
                buffer.differences[it->vertex()] -= 2.0*edge_costs[it->edge()];
            // edge to an element of the old set
            else if (lbl == old_label)
                buffer.differences[it->vertex()] += 2.0*edge_costs[it->edge()];
        }

        for (auto it = original_graph.adjacenciesFromVertexBegin(m.v); it != original_graph.adjacenciesFromVertexEnd(m.v); ++it)
        {
            if (buffer.is_moved[it->vertex()])
                continue;

            const auto lbl = buffer.vertex_labels[it->vertex()];

            // edge to an element of the new set
            if (lbl == m.new_label)
                --buffer.referenced_by[it->vertex()];
            // edge to an element of the old set
            else if (lbl == old_label)
            {
                ++buffer.referenced_by[it->vertex()];

                if (buffer.referenced_by[it->vertex()] == 1)
                    border.push_back(it->vertex());
            }
        }

        buffer.vertex_labels[m.v] = m.new_label;
        buffer.referenced_by[m.v] = 0;
        buffer.differences[m.v] = std::numeric_limits<double>::lowest();
        buffer.is_moved[m.v] = 1;
        moves.push_back(m);

        cumulative_diff += m.difference;

        if (cumulative_diff > max_move.first)
            max_move = std::make_pair(cumulative_diff, moves.size());
    }

    
    if (gain_from_merging > max_move.first && gain_from_merging > settings.epsilon)
    {
        A.insert(A.end(), B.begin(), B.end());

        for (auto a : A)
            buffer.vertex_labels[a] = label_A;

        for (auto b : B)
            buffer.vertex_labels[b] = label_A;

        B.clear();

        return gain_from_merging;
    }
    else if (max_move.first > settings.epsilon)
    {
        // revert some changes
        for (std::size_t i = max_move.second; i < moves.size(); ++i)
        {
            buffer.is_moved[moves[i].v] = 0;

            if (moves[i].new_label == label_B)
                buffer.vertex_labels[moves[i].v] = label_A;
            else
                buffer.vertex_labels[moves[i].v] = label_B;
        }

        // make sure that this is unique label
        if (B.empty())
            ++buffer.max_not_used_label;

        A.erase(std::partition(A.begin(), A.end(), [&](std::size_t a) { return !buffer.is_moved[a]; }), A.end());
        B.erase(std::partition(B.begin(), B.end(), [&](std::size_t b) { return !buffer.is_moved[b]; }), B.end());

        for (std::size_t i = 0; i < max_move.second; ++i)
            // move vertex to the other set
            if (moves[i].new_label == label_B)
                B.push_back(moves[i].v);
            else
                A.push_back(moves[i].v);

        return max_move.first;
    }
    else
        for (std::size_t i = 0; i < moves.size(); ++i)
            if (moves[i].new_label == label_B)
                buffer.vertex_labels[moves[i].v] = label_A;
            else
                buffer.vertex_labels[moves[i].v] = label_B;

    return .0;
}

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename SET, typename ECA>
inline double
kernighanLin(
    const ORIGINAL_GRAPH& original_graph,
    const LIFTED_GRAPH& lifted_graph,
    const ECA& edge_costs,
    SET& A,
    SET& B,
    const TwoCutSettings settings
) {
    TwoCutBuffers<ORIGINAL_GRAPH> buffer(original_graph);
    return kernighanLin(original_graph, lifted_graph, edge_costs, A, B, buffer, settings);
}

}
}
}

#endif
