#include <stdexcept>

#include "andres/ilp/gurobi.hxx"

#include "multicut.hxx"

int main() {
    typedef andres::ilp::Gurobi<> Ilp;

    testMulticut<Ilp>();

    return 0;
}
