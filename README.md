Graphs and Graph Algorithms in C++
============

Copyright (c) 2014 by Bjoern Andres, http://www.andres.sc/graph

This software was developed by Bjoern Andres, Duligur Ibeling, Giannis Kalofolias, Evgeny Levinkov and Mark Matten.
Enquiries shall be directed to bjoern@andres.sc.

Synopsis
------------

This set of header files implements directed and undirected graphs as adjacency lists with constant-time random access to neighboring vertices and incident edges.
Vertices and edges are always indexed by contiguous integers.
This indexing simplifies the implementation of algorithms for static graphs.
In dynamic settings where vertices and edges are removed from a graph,
indices of vertices and edges can change.
These changes can be followed, if necessary, via callbacks.
Subgraphs are defined by subgraph masks.


Features
------------

- Graphs and Digraphs, implemented as adjacency lists
- Constant-time random access to (neighboring) vertices and (incident) edges 
   - by contiguous integer indices
   - by STL-compliant random access iterators
- Insertion and removal of vertices and edges
- Multiple edges between the same pair of vertices (disabled by default)
- Callbacks for following changes of integer indices of vertices and edges
- Graph classes
   - Digraph
   - Graph
   - Complete graph
   - n-dimensional grid graph
- Graph algorithms
   - Connected components
      * by breadth-first search 
      * by disjoint sets
   - Shortest paths in weighted and unweighted graphs, as sequences of edges or vertices
      * Single source shortest path (SSSP)
      * Single pair shortest path (SPSP)
   - Cut-vertices and cut-edges (bridges)
   - Minimum spanning trees and minimum spanning forests
      * Prim's algorithm
      * Dynamic Programming
   - Maximum st-flow
      * Push-Relabel Algorithm with FIFO vertex selection rule
      * Edmonds-Karp Algorithm
   - Minimum Cost Multicut 
      * by integer programming, using Cplex or Gurobi
   - Set Partition
      * by a specialization of Minimum Cost Multicut for complete graphs


Developers
------------

- Bjoern Andres, bjoern@andres.sc
- Duligur Ibeling, duligur@gmail.com
- Giannis Kalofolias, kalofolias@mpi-inf.mpg.de
- Evgeny Levinkov, levinkov@mpi-inf.mpg.de
- Mark Matten, markmatten@gmail.com


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
