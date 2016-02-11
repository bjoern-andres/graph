#include <cstddef>

#include "andres/graph/graph.hxx"
#include "andres/graph/lifting.hxx"

#include <command-line-tools/probabilistic-lifting.hxx>

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

void testLiftEdgeValuesSPSP() {
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

    // define (two unit) edge values
    std::vector<unsigned int> edgeValues(graph.numberOfEdges(), 2);

    // lift graph
    const size_type distanceUpperBound = 4;
    andres::graph::Graph<> graphLifted;
    andres::graph::lift(graph, graphLifted, distanceUpperBound);

    // lift edge weights (by means of SPSP)
    std::vector<unsigned int> edgeValuesLifted(graphLifted.numberOfEdges());
    liftEdgeValues(graph, graphLifted, edgeValues.begin(), edgeValuesLifted.begin());

    // test lifted edge values
    for(size_type e = 0; e < graphLifted.numberOfEdges(); ++e)
    {
        const size_type v = graphLifted.vertexOfEdge(e, 0);
        const size_type w = graphLifted.vertexOfEdge(e, 1);
        const size_type distance = diff(v, w, graph.numberOfVertices());
        test(edgeValuesLifted[e] == distance * 2);
    }
}

void testLiftEdgeValuesFastMarching() {
    typedef std::size_t size_type;
    typedef andres::graph::GridGraph<2> GraphType;
    typedef GraphType::VertexCoordinate VertexCoordinate;
    typedef andres::graph::Graph<> LiftedGraphType;

    // build grid graph
    GraphType graph = {5, 4};

    // define (two unit) edge values
    std::vector<unsigned int> edgeValues(graph.numberOfEdges(), 2);

    // lift graph
    const size_type distanceUpperBound = 4;
    LiftedGraphType graphLifted;
    andres::graph::lift(graph, graphLifted, distanceUpperBound);

    // lift edge values (by menas of fast marching)
    const size_type interpolationOrder = 0;
    std::vector<unsigned int> edgeValuesLiftedFM(graphLifted.numberOfEdges());
    liftEdgeValues(graph, graphLifted, edgeValues.begin(), edgeValuesLiftedFM.begin(), interpolationOrder);

    // lift edge values (by menas of SPSP)
    std::vector<unsigned int> edgeValuesLiftedSPSP(graphLifted.numberOfEdges());
    liftEdgeValues<
        GraphType,
        LiftedGraphType,
        std::vector<unsigned int>::const_iterator,
        std::vector<unsigned int>::iterator
    >(
        graph,
        graphLifted,
        edgeValues.begin(),
        edgeValuesLiftedSPSP.begin()
    );

    // test lifted edge values
    VertexCoordinate cv, cw;
    for(size_type e = 0; e < graphLifted.numberOfEdges(); ++e)
        test(edgeValuesLiftedFM[e] == edgeValuesLiftedSPSP[e]);
}

int main() {
    testLiftGraph();
    testLiftGridGraphPathLengthMetric();
    testLiftGridGraphL2Metric();
    testLiftEdgeValuesSPSP();
    testLiftEdgeValuesFastMarching();

    return 0;
}
