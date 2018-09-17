#pragma once

#include <numeric>
#include "andres/partition.hxx"


namespace andres {
namespace graph {
namespace multicut {

namespace detail {

template<class MUTEX_STORAGE>
inline bool check_mutex(const std::size_t ru, const std::size_t rv,
                        const MUTEX_STORAGE & mutexes) {
    // get iterators to the mutex vectors of rep u and rep v
    auto mutex_it_u = mutexes[ru].begin();
    auto mutex_it_v = mutexes[rv].begin();

    // check if the mutex vectors contain the same mutex edge
    while (mutex_it_u != mutexes[ru].end() && mutex_it_v != mutexes[rv].end()) {
        if (*mutex_it_u < *mutex_it_v) {
            ++mutex_it_u;
        } else  {
            if (!(*mutex_it_v < *mutex_it_u)) {
                return true;
            }
            ++mutex_it_v;
        }
    }
    return false;
}


// insert 'dge_id' into the vectors containing mutexes of 'ru' and 'rv'
template<class MUTEX_STORAGE>
inline void insert_mutex(const std::size_t ru, const std::size_t rv,
                         const std::size_t edge_id, MUTEX_STORAGE & mutexes) {
    mutexes[ru].insert(std::upper_bound(mutexes[ru].begin(),
                                        mutexes[ru].end(),
                                        edge_id), edge_id);
    mutexes[rv].insert(std::upper_bound(mutexes[rv].begin(),
                                        mutexes[rv].end(),
                                        edge_id), edge_id);
}


// merge the mutex edges by merging from 'root_from' to 'root_to'
template<class MUTEX_STORAGE>
inline void merge_mutexes(const std::size_t root_from, const std::size_t root_to,
                          MUTEX_STORAGE & mutexes) {
    if (mutexes[root_from].size() == 0) {
        return;
    }

    if (mutexes[root_to].size() == 0){
        mutexes[root_to] = mutexes[root_from];
        return;
    }

    std::vector<std::size_t> merge_buffer;
    merge_buffer.reserve(std::max(mutexes[root_from].size(), mutexes[root_to].size()));

    std::merge(mutexes[root_from].begin(), mutexes[root_from].end(),
               mutexes[root_to].begin(), mutexes[root_to].end(),
               std::back_inserter(merge_buffer));

    mutexes[root_to] = merge_buffer;
    mutexes[root_from].clear();
}

}


// mutex watershed graph clustering
//
template<typename GRAPH, typename EVA, typename ELA>
void mutexWatershed(
    const GRAPH & graph,
    const EVA & edge_values,
    ELA & edge_labels
) {
    // sort the edge values in descending order
    // according to their absolute value
    const std::size_t n_edges = graph.numberOfEdges();
    const std::size_t n_vertices = graph.numberOfVertices();

    std::vector<std::size_t> indices(n_edges);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](const std::size_t a, const std::size_t b){
        return std::abs(edge_values[a]) > std::abs(edge_values[b]);
    });

    // data-structure storing persistent cuts a.k.a. mutexes
    typedef std::vector<std::vector<uint64_t>> MutexStorage;
    MutexStorage mutexes(n_vertices);

    // union find
    andres::Partition<std::size_t> partition(n_vertices);

    // iterate over the sorted edges
    for(const std::size_t edge_id : indices) {

        // check whether this edge is a repulsive edge (a.k.a mutex edge)
        // via the sign of its value
        // skip indifferent edges
        const auto edge_val = edge_values[edge_id];
        if(edge_val == 0)
            continue;
        const bool is_mutex_edge = edge_val < 0;

        // find the representatives of nodes connected by the edge
        const std::size_t u = graph.vertexOfEdge(edge_id, 0);
        const std::size_t v = graph.vertexOfEdge(edge_id, 1);
        std::size_t ru = partition.find(u);
        std::size_t rv = partition.find(v);

        // if the nodes are already connected, do nothing
        if(ru == rv) {
            continue;
        }

        // if we already have a mutex, we do not need to do anything
        // (if this is a regular edge, we do not link, if it is a mutex edge
        //  we do not need to insert the redundant mutex constraint)
        if(detail::check_mutex(ru, rv, mutexes)) {
            continue;
        }

        if(is_mutex_edge) {
            detail::insert_mutex(ru, rv, edge_id, mutexes);
        } else {
            // merge nodes and their mutex constraints
            partition.merge(u, v);
            // check  if we have to swap the roots
            if(partition.find(ru) == rv) {
                std::swap(ru, rv);
            }
            detail::merge_mutexes(rv, ru, mutexes);
        }
    }

    for (size_t i = 0; i < n_edges; ++i)
        edge_labels[i] = partition.find(graph.vertexOfEdge(i, 0)) == partition.find(graph.vertexOfEdge(i, 1)) ? 0 : 1;
}


}
}
}
