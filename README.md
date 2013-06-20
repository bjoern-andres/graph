Graphs and Graph Algorithms in C++
============

Copyright (c) 2013 by Bjoern Andres, http://www.andres.sc/


Synopsis
------------

This set of header files implements directed and undirected graphs as adjacency lists.
Vertices and edges are always indexed by contiguous integers.
This indexing simplifies the implementation of algorithms for static graphs.
In dynamic settings where vertices and edges are removed from a graph,
indices of vertices and edges can change.
These changes can be followed, if necessary, by means of a visitor.
Subgraphs are defined by subgraph masks.


Features
------------

- Directed and undirected graphs implemented as adjacency lists with constant-time access
- Access to vertices and edges by contiguous integer indices
- Access to neighboring vertices and incident edges by STL-compliant random access iterators
- Insertion and removal of vertices and edges
- Multiple edges, which are disabled by default, can be enabled
- Visitors that follow changes of vertex and edge indices
- Algorithms
  - Connected components by breadth-first search and disjoint sets
  - Shortest paths (SSSP, SPSP) in weighted and unweighted graphs
  - Maximum s-t-flow (Push-Relabel Algorithm with FIFO vertex selection rule)
  - Minimum multicuts by interger programming, using Cplex or Gurobi


Contributions
------------

- Mark Matten, markmatten@gmail.com
  - Push-Relabel Algorithm with FIFO vertex selection rule for computing maximum s-t-flow


License
------------

Copyright (c) 2013 by Bjoern Andres, http://www.andres.sc/

This software was developed by Bjoern Andres.
Enquiries shall be directed to bjoern@andres.sc.

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
