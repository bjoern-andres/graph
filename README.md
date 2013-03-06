Graph
============

Graphs with Integer Indexing and Graph Algorithms in C++


Short description

This set of header files implements directed and undirected graphs as adjacency lists.
Vertices and edges are always indexed by contiguous integers.
This indexing simplifies the implementation of algorithms for static graphs.
In dynamic settings where vertices and edges are removed from a graph,
indices of vertices and edges can change.
These changes can be followed, if necessary, by means of a visitor.
Subgraphs are defined by subgraph masks.


Features

- Directed and undirected graphs, implemented as adjacency lists
- Access to vertices and edges by contiguous integer indices
- Access to neighboring vertices and incident edges by STL-compliant random access iterators
- Insertion and removal of vertices and edges
- Multiple edges, which are disabled by default, can be enabled
- Visitors that follow changes of vertex and edge indices
- Algorithms
  - Minimum multicuts by interger programming, using Cplex or Gurobi
  - Connected components by breadth-first search and disjoint sets
  - Shortest paths (SSSP, SPSP) in weighted and unweighted graphs
