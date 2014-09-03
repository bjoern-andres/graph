#include <stdexcept>

#include "andres/ilp/cplex.hxx"

#include "multicut.hxx"

int main() {
    typedef andres::ilp::Cplex<> Ilp;

    testMulticut<Ilp>();
    testMulticutCompleteGraph<Ilp>();

    return 0;
}
