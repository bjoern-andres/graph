#pragma once
#ifndef ANDRES_GRAPH_DETAIL_BRIDGES_HXX
#define ANDRES_GRAPH_DETAIL_BRIDGES_HXX

namespace andres {
namespace graph {
namespace detail {

template<typename GRAPH, typename SUBGRAPH, typename BUFFER>
inline
void find_bridges(const GRAPH& graph, const SUBGRAPH& subgraph_mask, std::size_t starting_vertex, BUFFER& buffer)
{
    std::fill(buffer.visited_.begin(), buffer.visited_.end(), 0);

    std::stack<std::size_t> S;

    S.push(starting_vertex);
    buffer.depth_[starting_vertex] = 0;
    buffer.parent_[starting_vertex] = -1;

    while (!S.empty())
    {
        auto v = S.top();
        S.pop();

        if (!buffer.visited_[v])
        {
            buffer.visited_[v] = 1;
            buffer.next_out_arc_[v] = graph.verticesFromVertexBegin(v);
            buffer.min_successor_depth_[v] = buffer.depth_[v];
        }
        else
        {
            auto to = *buffer.next_out_arc_[v];

            if (buffer.min_successor_depth_[to] > buffer.depth_[v])
            {
                typename GRAPH::EdgeIterator e = graph.edgesFromVertexBegin(v) + (buffer.next_out_arc_[v] - graph.verticesFromVertexBegin(v));
                buffer.is_bridge_[*e] = 1;
            }

            buffer.min_successor_depth_[v] = std::min(buffer.min_successor_depth_[v], buffer.min_successor_depth_[to]);
            ++buffer.next_out_arc_[v];
        }

        while (buffer.next_out_arc_[v] != graph.verticesFromVertexEnd(v))
        {
            typename GRAPH::EdgeIterator e = graph.edgesFromVertexBegin(v) + (buffer.next_out_arc_[v] - graph.verticesFromVertexBegin(v));

            if (
                !subgraph_mask.vertex(*buffer.next_out_arc_[v]) ||
                !subgraph_mask.edge(*e)
                )
            {
                ++buffer.next_out_arc_[v];
                continue;
            }

            auto to = *buffer.next_out_arc_[v];
            if (buffer.visited_[to])
            {
                if(buffer.parent_[v] != to)
                    buffer.min_successor_depth_[v] = std::min(buffer.min_successor_depth_[v], buffer.depth_[to]);

                ++buffer.next_out_arc_[v];
            }
            else
            {
                S.push(v);
                S.push(to);
                buffer.parent_[to] = v;
                buffer.depth_[to] = buffer.depth_[v] + 1;
                buffer.is_bridge_[*e] = 0;
                break;
            }
        }
    }
}

}
}
}
#endif
