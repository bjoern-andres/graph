/// Copyright (c) 2013 by Mark Matten
/// 
/// This software was developed by Mark Matten.
/// Enquiries shall be directed to markmatten@gmail.com
///
/// Redistribution and use in source and binary forms, with or without 
/// modification, are permitted provided that the following conditions are met:
///
/// - Redistributions of source code must retain the above copyright notice,
///   this list of conditions and the following disclaimer.
/// - Redistributions in binary form must reproduce the above copyright notice, 
///   this list of conditions and the following disclaimer in the documentation
///   and/or other materials provided with the distribution.
/// - The name of the author must not be used to endorse or promote products 
///   derived from this software without specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
/// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
/// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
/// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
/// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
/// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
/// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
/// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
/// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///
#pragma once
#ifndef ANDRES_GRAPH_MAX_FLOW_HXX
#define ANDRES_GRAPH_MAX_FLOW_HXX

#include <cassert>
#include <vector>
#include <queue>
#include <limits> // std::numeric_limit
#include <algorithm> // std::max, std::min

namespace andres {
namespace graph {

template<class GRAPH, class FLOW>
class MaxFlowPushRelabel {
public:
	typedef GRAPH GraphType;
	typedef FLOW Flow;
	typedef typename GraphType::VertexIterator VertexIterator;
	typedef typename GraphType::EdgeIterator EdgeIterator;
	typedef typename GraphType::AdjacencyIterator AdjacencyIterator;

	MaxFlowPushRelabel();
	void clear();
	Flow maxFlow() const;
	Flow flow(const size_t) const;
	size_t numberOfPushes() const;
	size_t numberOfRelabels() const;
	template<class EDGE_WEIGHT_ITERATOR>
		MaxFlowPushRelabel(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);
	template<class EDGE_WEIGHT_ITERATOR>
		Flow operator()(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t, const size_t);

private:
	template<class EDGE_WEIGHT_ITERATOR>
		void push(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t);
	template<class EDGE_WEIGHT_ITERATOR>
		void pushBack(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t);
	template<class EDGE_WEIGHT_ITERATOR>
		void relabel(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t);
	template<class EDGE_WEIGHT_ITERATOR>
		void discharge(const GraphType&, EDGE_WEIGHT_ITERATOR, const size_t);
	void gapRelabel(const GraphType&, const size_t);

	std::vector<size_t> height_;
	std::vector<size_t> labelCount_;	
	std::vector<Flow> excess_;
	std::vector<Flow> flow_;
	std::vector<bool> active_;
	std::queue<size_t> queue_;
	size_t sourceVertexIndex_;
	size_t sinkVertexIndex_;
	size_t pushCount_;
	size_t relabelCount_;
};

template<class GRAPH, class FLOW>
inline
MaxFlowPushRelabel<GRAPH, FLOW>::MaxFlowPushRelabel()
:	height_(), 
	labelCount_(),
	excess_(), 
	flow_(),
	active_(),
	queue_(),
	sourceVertexIndex_(),
	sinkVertexIndex_(),
	pushCount_(), 
	relabelCount_()
{}

template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline
MaxFlowPushRelabel<GRAPH, FLOW>::MaxFlowPushRelabel(
	const GraphType& graph,
	EDGE_WEIGHT_ITERATOR edgeWeightIterator,
	const size_t sourceVertexIndex,
	const size_t sinkVertexIndex
)
:	height_(), 
	labelCount_(),
	excess_(), 
	flow_(),
	active_(),
	queue_(),
	sourceVertexIndex_(),
	sinkVertexIndex_(),
	pushCount_(), 
	relabelCount_()
{
	(*this)(graph, edgeWeightIterator, sourceVertexIndex, sinkVertexIndex);
}

template<class GRAPH, class FLOW>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::clear() {
	height_.clear();
	labelCount_.clear();
	excess_.clear();
	flow_.clear();
	active_.clear();
	sourceVertexIndex_ = 0;
	sinkVertexIndex_ = 0;
	pushCount_ = 0;
	relabelCount_ = 0;

	assert(queue_.empty());
}

template<class GRAPH, class FLOW>
inline typename MaxFlowPushRelabel<GRAPH, FLOW>::Flow 
MaxFlowPushRelabel<GRAPH, FLOW>::maxFlow() const  {
	return excess_[sinkVertexIndex_];
}

template<class GRAPH, class FLOW>
inline typename MaxFlowPushRelabel<GRAPH, FLOW>::Flow 
MaxFlowPushRelabel<GRAPH, FLOW>::flow(
	const size_t edgeIndex
) const {
	assert(edgeIndex < flow_.size());
	return flow_[edgeIndex];
}

template<class GRAPH, class FLOW>
inline size_t 
MaxFlowPushRelabel<GRAPH, FLOW>::numberOfPushes() const {
	return pushCount_;
}

template<class GRAPH, class FLOW>
inline size_t 
MaxFlowPushRelabel<GRAPH, FLOW>::numberOfRelabels() const {
	return relabelCount_;
}

template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline typename MaxFlowPushRelabel<GRAPH, FLOW>::Flow 
MaxFlowPushRelabel<GRAPH, FLOW>::operator()(
	const GraphType& graph,
	EDGE_WEIGHT_ITERATOR edgeWeightIterator,
	const size_t sourceVertexIndex,
	const size_t sinkVertexIndex
) {
	const size_t numberOfVertices = graph.numberOfVertices();
	const size_t numberOfEdges = graph.numberOfEdges();

	// initialization
	sourceVertexIndex_ = sourceVertexIndex;
	sinkVertexIndex_ = sinkVertexIndex;
	pushCount_ = 0;
	relabelCount_ = 0;
	height_.resize(numberOfVertices);
	excess_.resize(numberOfVertices);
	active_.resize(numberOfVertices);
	std::fill(height_.begin(), height_.end(), size_t());
	std::fill(excess_.begin(), excess_.end(), Flow());
	std::fill(active_.begin(), active_.end(), false);
	height_[sourceVertexIndex] = numberOfVertices;
	excess_[sourceVertexIndex] = std::numeric_limits<size_t>::max(); // this is supposed to be infinite
	active_[sourceVertexIndex] = true;
	active_[sinkVertexIndex] = true;
	labelCount_.resize((2 * numberOfVertices) + 2); // 2n + 1 is the maximum possible height for a vertex
	std::fill(labelCount_.begin(), labelCount_.end(), size_t());
	labelCount_[0] = numberOfVertices - 1;
	labelCount_[numberOfVertices] = 1;
	flow_.resize(numberOfEdges);
	std::fill(flow_.begin(), flow_.end(), Flow());

	// first, push as much flow as possible from the source to all adjacent vertices
	for(EdgeIterator it = graph.edgesFromVertexBegin(sourceVertexIndex); it != graph.edgesFromVertexEnd(sourceVertexIndex); ++it) {
		const size_t edgeIndex = *it;
		push(graph, edgeWeightIterator, edgeIndex);
	}
	
	while(!queue_.empty()) { // main loop
		const size_t u = queue_.front();
		queue_.pop();
		active_[u] = false;
		discharge(graph, edgeWeightIterator, u);
		active_[u] = false; // TODO: why does putting active_[u] = false twice decrease the number of pushes?
	}

	return maxFlow();
}

template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::push(
	const GRAPH& graph,
	EDGE_WEIGHT_ITERATOR edgeWeightIterator,
	const size_t edgeIndex
) {
	const size_t u = graph.vertexOfEdge(edgeIndex, 0);
	const size_t v = graph.vertexOfEdge(edgeIndex, 1);
	// push from u to v
	Flow pushAmount = std::min(excess_[u], edgeWeightIterator[edgeIndex] - flow_[edgeIndex]);
	flow_[edgeIndex] += pushAmount;
	excess_[u] -= pushAmount;
	excess_[v] += pushAmount;
	if (!active_[v] && excess_[v] > 0) {
		active_[v] = true;
		queue_.push(v);
	}
	pushCount_++;
}

template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::pushBack(
	const GRAPH& graph,
	EDGE_WEIGHT_ITERATOR edgeWeightIterator,
	const size_t edgeIndex
) {
	const size_t u = graph.vertexOfEdge(edgeIndex, 0);
	const size_t v = graph.vertexOfEdge(edgeIndex, 1);
	// push back from v to u
	Flow pushAmount = std::min(excess_[v], flow_[edgeIndex]);
	flow_[edgeIndex] -= pushAmount;
	excess_[u] += pushAmount;
	excess_[v] -= pushAmount;
	if (!active_[u] && excess_[u] > 0) {
		active_[u] = true;
		queue_.push(u);
	}
	pushCount_++;
}

/// increase height of u to 1 greater than the minimum height of its neighbors.
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::relabel(
	const GRAPH& graph,
	EDGE_WEIGHT_ITERATOR edgeWeightIterator,
	const size_t u
) {
	size_t minHeight = (2 * graph.numberOfVertices()) + 1;
	const size_t oldHeight = height_[u];
	for(EdgeIterator it = graph.edgesFromVertexBegin(u); it != graph.edgesFromVertexEnd(u); ++it) {
		const size_t edgeIndex = *it;
		const size_t v = graph.vertexOfEdge(edgeIndex, 1); // edge is (u, v)
		if(edgeWeightIterator[edgeIndex] - flow_[edgeIndex] > 0) {
			minHeight = std::min(minHeight, height_[v]);
		}
	}
	for(EdgeIterator it = graph.edgesToVertexBegin(u); it != graph.edgesToVertexEnd(u); ++it) {
		const size_t edgeIndex = *it;
		const size_t v = graph.vertexOfEdge(edgeIndex, 0); // edge is (v, u)
		if(flow_[edgeIndex] > 0) {
			minHeight = std::min(minHeight, height_[v]);
		}
	}
	height_[u] = minHeight + 1;
	if (!active_[u] && excess_[u] > 0) {
		active_[u] = true;
		queue_.push(u);
	}
	relabelCount_++;

	// gap relabeling
	labelCount_[oldHeight]--;
	labelCount_[height_[u]]++;
	if (labelCount_[oldHeight] == 0) {
		gapRelabel(graph, oldHeight);
	}
}

/// while there is excess flow at u, try to push flow to neighbors. If no push is available, relabel u.
template<class GRAPH, class FLOW>
template<class EDGE_WEIGHT_ITERATOR>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::discharge(
	const GRAPH& graph,
	EDGE_WEIGHT_ITERATOR edgeWeightIterator,
	const size_t u
) {
	size_t seenForward = 0, seenBackward = 0;
	while(excess_[u] > 0) {
		for(EdgeIterator it = graph.edgesFromVertexBegin(u); it != graph.edgesFromVertexEnd(u); ++it) {
			const size_t edgeIndex = *it;
			const size_t v = graph.vertexOfEdge(edgeIndex, 1); // edge is the pair (u, v)
			if(edgeWeightIterator[edgeIndex] - flow_[edgeIndex] > 0 && height_[u] > height_[v]) {
				push(graph, edgeWeightIterator, edgeIndex);
			}
		}
		for(EdgeIterator it = graph.edgesToVertexBegin(u); it != graph.edgesToVertexEnd(u); ++it) {
			const size_t edgeIndex = *it;
			const size_t v = graph.vertexOfEdge(edgeIndex, 0); // edge is the pair (v, u)
			if(flow_[edgeIndex] > 0 && height_[u] > height_[v]) {
				pushBack(graph, edgeWeightIterator, edgeIndex);
			}
		}
		relabel(graph, edgeWeightIterator, u);
	}
}

template<class GRAPH, class FLOW>
inline void
MaxFlowPushRelabel<GRAPH, FLOW>::gapRelabel(
	const GRAPH& graph,
	const size_t threshold
) {
	for(size_t i = 0; i < graph.numberOfVertices(); i++) {
		if(height_[i] > threshold) {
			height_[i] = std::max(height_[i], graph.numberOfVertices() + 1);;
		}
	}
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MAX_FLOW_HXX
