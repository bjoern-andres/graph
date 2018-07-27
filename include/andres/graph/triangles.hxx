#pragma once
#ifndef ANDRES_GRAPH_TRIANGLES_HXX
#define ANDRES_GRAPH_TRIANGLES_HXX

#include <list>

namespace andres {
namespace graph {

// Efficiently list all triangles of a graph.
// Output: Vector of vertex triples    
//
// Implementation of "forward" algorithm from [1]
//
// [1] Schank, Wagner. Finding, Counting and Listing all
// Triangles in Large Graphs, An Experimental Study. 2005
//
template<typename GRAPH>
std::vector<std::vector<size_t>> findTriangles(const GRAPH& graph)
{
    std::vector<std::vector<size_t>> triangles;

    // get degrees of all vertices
    std::vector<size_t> vertices(graph.numberOfVertices());
    std::vector<size_t> index(graph.numberOfVertices());
    std::vector<size_t> degrees(graph.numberOfVertices());
    for (size_t v = 0; v < graph.numberOfVertices(); v++)
    {
        vertices[v] = v;
        degrees[v] = graph.numberOfEdgesFromVertex(v);
    }

    // sort vertices according to degrees
    std::sort(vertices.begin(), vertices.end(), [&] (const size_t& a, const size_t& b) {return (degrees[a] < degrees[b]);});

    // determine mapping of vertices to indices
    for (size_t i = 0; i < graph.numberOfVertices(); i++)
    {
        auto v = vertices[i];
        index[v] = i;
    }

    // dynamic adjacency data
    std::vector<std::list<size_t>> adj(graph.numberOfVertices());

    // MAIN LOOP
    for (size_t i = 0; i < graph.numberOfVertices(); i++)
    {
        auto u = vertices[i];
        // iterate over neighbors
        for (auto it = graph.adjacenciesFromVertexBegin(u); it != graph.adjacenciesFromVertexEnd(u); it++)
        {
            auto v = it->vertex();
            auto j = index[v];

            // only consider smaller indices
            if (i < j)
            {
                // find intersection in adjacency lists
                auto it_u = adj[i].begin();
                auto it_v = adj[j].begin();

                while (it_u != adj[i].end() && it_v != adj[j].end())
                {
                    if (*it_u == *it_v)
                    {
                        auto k = *it_u;
                        std::vector<size_t> triangle{u,v,vertices[k]}; 
                        triangles.push_back(triangle);
                        it_u++;
                        it_v++;
                    }
                    else if (*it_u < *it_v)
                        it_u++;
                    else if (*it_v < *it_u)
                        it_v++;
                }

                adj[j].push_back(i);
            }
        }
    }

    return triangles;
}
    
}
}

#endif
