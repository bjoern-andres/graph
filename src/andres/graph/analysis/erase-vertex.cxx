#include <random>
#include <vector>
#include <chrono>
#include <iostream>

#include <andres/graph/graph.hxx>

int main() {
    // construct graph
    size_t const numberOfVertices = 1000;
    andres::graph::Graph<> graph(numberOfVertices);
    for(size_t j = 0; j < numberOfVertices; ++j)
    for(size_t k = 0; k < j; ++k) {
        graph.insertEdge(j, k);
    }

    // remove verices
    typedef std::chrono::high_resolution_clock Clock;
    std::default_random_engine randomEngine;
    auto const timeStart = Clock::now();
    for(size_t j = 0; j < numberOfVertices; ++j) {
        std::uniform_int_distribution<size_t> distribution(0, graph.numberOfVertices() - 1);
        size_t const index = distribution(randomEngine);
        graph.eraseVertex(index);
    }
    auto const timeStop = Clock::now();
    auto const duration = std::chrono::duration<double, std::milli>(timeStop - timeStart);
    std::cout << "erasing " << numberOfVertices
              << " vertices has taken " << duration.count()
              << " ms." << std::endl;

    return 0;
}
