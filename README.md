Graphs and Graph Algorithms in C++
============

Copyright (c) 2014 by Bjoern Andres, bjoern@andres.sc

http://www.andres.sc/graph.html


Developers
------------

- Bjoern Andres, bjoern@andres.sc
- Duligur Ibeling, duligur@gmail.com
- Giannis Kalofolias, kalofolias@mpi-inf.mpg.de
- Evgeny Levinkov, levinkov@mpi-inf.mpg.de
- Mark Matten, markmatten@gmail.com


Synopsis
------------

This set of header files implements directed and undirected graphs as adjacency lists with constant-time random access to neighboring vertices and incident edges.
Vertices and edges are always indexed by contiguous integers.
This indexing simplifies the implementation of algorithms for static graphs.
In dynamic settings where vertices and edges are removed from a graph,
indices of vertices and edges can change.
These changes can be followed, if necessary, via a functor.
Subgraphs are defined by subgraph masks.


Features
------------

- Graphs and Digraphs, implemented as adjacency lists
   - Subgraphs defined by subgraph masks
   - Multiple edges between the same pair of vertices (disabled by default)
   - Insertion and removal of vertices and edges
   - Tracking of integer indices of vertices and edges via a functor 
   - Constant-time random access to neighboring vertices and incident edges 
      - By contiguous integer indices
      - By STL-compliant random access iterators
- Graph classes
   - Digraph
   - Graph
   - Complete graph
   - n-dimensional grid graph
- Graph file format and I/O based on HDF5
- Graph algorithms
   - Search and traversal with functors
      * depth first search (DFS)
      * breadth first search (BFS)
   - 1-connected components
      * by breadth first search 
      * by disjoint sets
   - 2-connected components
      * cut-vertices/articulation vertices (Hopcroft-Tarjan Algorithm)
      * cut-edges/bridges (Tarjan's algorithm)
   - Shortest paths in weighted and unweighted graphs, as sequences of edges or vertices
      * Single source shortest path (SSSP)
      * Single pair shortest path (SPSP)
   - Minimum spanning trees and minimum spanning forests
      * Prim's algorithm
      * Dynamic Programming
   - Maximum st-flow
      * Goldberg-Tarjan Algorithm/Push-Relabel Algorithm with FIFO vertex selection rule
      * Edmonds-Karp Algorithm
   - Minimum Cost Multicut 
      * by integer programming, using Cplex or Gurobi
   - Set Partition
      * by a specialization of Minimum Cost Multicut for complete graphs


References
------------

- S. Chopra and M. R. Rao. The partition problem. Mathematical Programming 59(1-3):87-115. 1993
- A. V. Goldberg and R. E. Tarjan. A new approach to the maximum-flow problem. Journal of the ACM 35(4):921-940. 1988
- R. E. Tarjan. A note on finding the bridges of a graph. Information Processing Letters 2(6):160-161. 1974
- J. Hopcroft and R. E. Tarjan. Efficient algorithms for graph manipulation. Communications of the ACM 16(6):372-378. 1973
- J. Edmonds and R. M. Karp. Theoretical improvements in algorithmic efficiency for network flow problems. Journal of the ACM 19(2):248-264. 1972
- R. C. Prim. Shortest connection networks and some generalizations. Bell System Technical Journal 36:1389-1401. 1957


License
------------

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, 
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- The name of the author must not be used to endorse or promote products 
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
