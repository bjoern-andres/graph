#include <cstddef>

#include "andres/graph/graph.hxx"
#include "andres/graph/lifting.hxx"


inline void test(const bool condition) {
    if(!condition) throw std::logic_error("test failed.");
}

// difference modulo m
template<class T>
inline T diff(const T a, const T b, const T m) {
    const T d = a < b ? b - a : a - b;
    const T d2 = m - d;
    return d < d2 ? d : d2;
}

void testLiftGraph() {
    typedef std::size_t size_type;

    // build cycle
    const size_type numberOfVertices = 10;
    andres::graph::Graph<> graph;
    graph.insertVertices(numberOfVertices);
    for(size_type v = 0; v < graph.numberOfVertices(); ++v)
    {
        const int w = (v + 1) % graph.numberOfVertices();
        graph.insertEdge(v, w);
    }

    for(size_type distanceLowerBound = 0; distanceLowerBound < graph.numberOfVertices(); ++distanceLowerBound)
        for(size_type distanceUpperBound = 0; distanceUpperBound < graph.numberOfVertices(); ++distanceUpperBound)
        {
            andres::graph::Graph<> graphLifted;
            andres::graph::lift(graph, graphLifted, distanceUpperBound, distanceLowerBound);
            test(graphLifted.numberOfVertices() == graph.numberOfVertices());
            for(size_type v = 0; v < graphLifted.numberOfVertices(); ++v)
                for(size_type w = 0; w < graphLifted.numberOfVertices(); ++w)
                {
                    std::pair<bool, std::size_t> p = graphLifted.findEdge(v, w);
                    const size_type d = diff(v, w, graph.numberOfVertices());
                    const bool withinRange = (d > distanceLowerBound && d <= distanceUpperBound);
                    test(p.first == withinRange);
                }
        }
}

void testLiftGridGraphPathLengthMetric() {
    typedef std::size_t size_type;

    andres::graph::GridGraph<2> gridGraph = {5, 4};
    andres::graph::Graph<> graph(gridGraph.numberOfVertices());
    graph.reserveEdges(gridGraph.numberOfEdges());
    for(size_type e = 0; e < gridGraph.numberOfEdges(); ++e)
        graph.insertEdge(gridGraph.vertexOfEdge(e, 0), gridGraph.vertexOfEdge(e, 1));

    for(size_type distanceLowerBound = 0; distanceLowerBound < graph.numberOfVertices(); ++distanceLowerBound)
        for(size_type distanceUpperBound = 0; distanceUpperBound < graph.numberOfVertices(); ++distanceUpperBound)
        {
            andres::graph::Graph<> gridGraphLifted;
            andres::graph::lift(gridGraph, gridGraphLifted, distanceUpperBound, distanceLowerBound, andres::graph::LiftingMetric::PathLength);

            andres::graph::Graph<> graphLifted;
            andres::graph::lift(graph, graphLifted, distanceUpperBound, distanceLowerBound); // tested by function testLiftGraph

            test(gridGraphLifted.numberOfVertices() == gridGraph.numberOfVertices());
            test(gridGraphLifted.numberOfEdges() == graphLifted.numberOfEdges());
            for(size_type v = 0; v < graphLifted.numberOfVertices(); ++v)
                for(size_type w = 0; w < graphLifted.numberOfVertices(); ++w)
                {
                    std::pair<bool, std::size_t> p1 = gridGraphLifted.findEdge(v, w);
                    std::pair<bool, std::size_t> p2 = graphLifted.findEdge(v, w);
                    test(p1.first == p2.first);
                }
        }
}

void testLiftGridGraphL2Metric() {
    typedef std::size_t size_type;
    typedef andres::graph::GridGraph<2> GraphType;
    typedef GraphType::VertexCoordinate VertexCoordinate;
    typedef andres::graph::Graph<> LiftedGraphType;

    GraphType graph = {5, 4};
    for(size_type distanceLowerBound = 0; distanceLowerBound < graph.numberOfVertices(); ++distanceLowerBound)
        for(size_type distanceUpperBound = 0; distanceUpperBound < graph.numberOfVertices(); ++distanceUpperBound)
        {
            const size_type distanceLowerBoundSquared = distanceLowerBound * distanceLowerBound;
            const size_type distanceUpperBoundSquared = distanceUpperBound * distanceUpperBound;

            LiftedGraphType graphLifted;
            andres::graph::lift(graph, graphLifted, distanceUpperBound, distanceLowerBound, andres::graph::LiftingMetric::L2);
            /*
            std::cout << distanceLowerBound << " < d < " << distanceUpperBound << std::endl;
            for(size_type e = 0; e < graphLifted.numberOfEdges(); ++e) {
                std::cout << graphLifted.vertexOfEdge(e, 0) << " " << graphLifted.vertexOfEdge(e, 1) << std::endl;
            }
            */
            test(graphLifted.numberOfVertices() == graph.numberOfVertices());
            for(size_type v = 0; v < graphLifted.numberOfVertices(); ++v)
            {
                VertexCoordinate cv;
                graph.vertex(v, cv);

                for(size_type w = 0; w < graphLifted.numberOfVertices(); ++w)
                {
                    VertexCoordinate cw;
                    graph.vertex(w, cw);

                    size_type distanceSquared = 0;
                    for(size_type j = 0; j < GraphType::DIMENSION; ++j)
                    {
                        const size_type d = cv[j] - cw[j];
                        distanceSquared += d * d;
                    }

                    const bool withinRange = (distanceSquared > distanceLowerBoundSquared && distanceSquared <= distanceUpperBoundSquared);

                    std::pair<bool, std::size_t> p = graphLifted.findEdge(v, w);
                    test(p.first == withinRange);
                }
            }
        }
}

int main() {
    testLiftGraph();
    testLiftGridGraphPathLengthMetric();
    testLiftGridGraphL2Metric();

    return 0;
}
