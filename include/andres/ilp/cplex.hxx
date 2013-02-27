// Copyright (c) 2013 by Bjoern Andres.
// 
// This software was developed by Bjoern Andres.
// Enquiries shall be directed to bjoern@andres.sc.
//
// All advertising materials mentioning features or use of this software must
// display the following acknowledgement: ``This product includes andres::ilp
// developed by Bjoern Andres. Please direct enquiries concerning andres::ilp
// to bjoern@andres.sc''.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// - All advertising materials mentioning features or use of this software must 
//   display the following acknowledgement: ``This product includes 
//   andres::ilp developed by Bjoern Andres. Please direct enquiries 
//   concerning andres::ilp to bjoern@andres.sc''.
// - The name of the author must not be used to endorse or promote products 
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#pragma once
#ifndef ANDRES_CPLEX_HXX
#define ANDRES_CPLEX_HXX

#include "ilcplex/ilocplex.h"

namespace andres {
namespace ilp {

template<class T = double>
class Cplex {
public:
    typedef T value_type;
    enum PreSolver {PRE_SOLVER_AUTO, PRE_SOLVER_PRIMAL, PRE_SOLVER_DUAL, PRE_SOLVER_NONE};
    enum LPSolver {LP_SOLVER_AUTO, LP_SOLVER_PRIMAL_SIMPLEX, LP_SOLVER_DUAL_SIMPLEX, LP_SOLVER_NETWORK_SIMPLEX, LP_SOLVER_BARRIER, LP_SOLVER_SIFTING, LP_SOLVER_CONCURRENT};

    Cplex();
    ~Cplex();
    void setNumberOfThreads(const size_t);
    void setAbsoluteGap(const value_type);
    void setRelativeGap(const value_type);
    void setVerbosity(const bool);
    void initModel(const size_t, const value_type*);
    void setLPSolver(const LPSolver);
    void setPreSolver(const PreSolver, const int = -1);
    template<class VariableIndexIterator, class CoefficientIterator>
        void addConstraint(VariableIndexIterator, VariableIndexIterator,
                           CoefficientIterator, const value_type, const value_type);
    template<class Iterator>
        void setStart(Iterator);

    void optimize();

    value_type label(const size_t) const;
    size_t numberOfThreads() const;
    value_type absoluteGap() const;
    value_type relativeGap() const;

private:
    IloEnv ilpEnvironment_;
    IloObjective ilpObjective_;
    IloModel ilpModel_;
    IloCplex ilpSolver_;
    IloNumArray ilpSolution_;
    IloNumVarArray ilpVariables_;
    IloNumArray ilpObjectiveCoefficients_;
    IloNumArray ilpStartValues_;
};

template<class T>
inline
Cplex<T>::Cplex()
:   ilpEnvironment_(),
    ilpObjective_(ilpEnvironment_),
    ilpModel_(ilpEnvironment_),
    ilpSolver_(ilpModel_),
    ilpSolution_(ilpEnvironment_),
    ilpVariables_(),
    ilpObjectiveCoefficients_(),
    ilpStartValues_()
{
    setVerbosity(false);
}

template<class T>
Cplex<T>::~Cplex() {
    ilpEnvironment_.end();
}

template<class T>
inline void
Cplex<T>::setNumberOfThreads (
    const size_t numberOfThreads
) {
    ilpSolver_.setParam(IloCplex::Threads, numberOfThreads);
}

template<class T>
inline void
Cplex<T>::setAbsoluteGap(
    const T gap
) {
    ilpSolver_.setParam(IloCplex::EpAGap, gap);
}

template<class T>
inline void
Cplex<T>::setRelativeGap(
    const T gap
) {
    ilpSolver_.setParam(IloCplex::EpGap, gap);
}

template<class T>
inline void
Cplex<T>::setVerbosity(
    const bool verbosity
) {
    if(verbosity) {
        ilpSolver_.setParam(IloCplex::MIPDisplay, 2);
        ilpSolver_.setParam(IloCplex::BarDisplay, 1);
        ilpSolver_.setParam(IloCplex::NetDisplay, 1);
        ilpSolver_.setParam(IloCplex::SiftDisplay, 1);
        ilpSolver_.setParam(IloCplex::SimDisplay, 1);
    }
    else {
        ilpSolver_.setParam(IloCplex::MIPDisplay, 0);
        ilpSolver_.setParam(IloCplex::BarDisplay, 0);
        ilpSolver_.setParam(IloCplex::NetDisplay, 0);
        ilpSolver_.setParam(IloCplex::SiftDisplay, 0);
        ilpSolver_.setParam(IloCplex::SimDisplay, 0);
    }
}

template<class T>
inline void
Cplex<T>::initModel(
    const size_t numberOfVariables,
    const T* coefficients
) {
    ilpObjective_.setSense(IloObjective::Minimize);
    ilpVariables_ = IloNumVarArray(ilpEnvironment_, numberOfVariables, 0, 1, ILOBOOL);
    ilpObjectiveCoefficients_ = IloNumArray(ilpEnvironment_, numberOfVariables);
    for(size_t j=0; j<numberOfVariables; ++j) {
        ilpObjectiveCoefficients_[j] = coefficients[j];
    }
    ilpObjective_.setLinearCoefs(ilpVariables_, ilpObjectiveCoefficients_);
    ilpModel_.add(ilpObjective_);
    ilpStartValues_ = IloNumArray(ilpEnvironment_, numberOfVariables);
    for(size_t j=0; j<numberOfVariables; ++j) {
        ilpStartValues_[j] = 1;
    }
}

template<class T>
inline void
Cplex<T>::optimize() {
    ilpSolver_.solve();
    ilpSolver_.getValues(ilpSolution_, ilpVariables_);
}

template<class T>
inline T
Cplex<T>::label (
    const size_t variableIndex
) const {
    return ilpSolution_[variableIndex];
}

template<class T>
inline size_t
Cplex<T>::numberOfThreads() const {
    return ilpSolver_.getParam(IloCplex::Threads);
}

template<class T>
inline T
Cplex<T>::absoluteGap() const {
    return ilpSolver_.getParam(IloCplex::EpAGap);
}

template<class T>
inline T
Cplex<T>::relativeGap() const {
    return ilpSolver_.getParam(IloCplex::EpGap);
}

template<class T>
template<class VariableIndexIterator, class CoefficientIterator>
void Cplex<T>::addConstraint(
    VariableIndexIterator viBegin,
    VariableIndexIterator viEnd,
    CoefficientIterator coefficient,
    const T lowerBound,
    const T upperBound
) {
    IloRange constraint(ilpEnvironment_, lowerBound, upperBound);
    for(; viBegin != viEnd; ++viBegin, ++coefficient) {
       constraint.setLinearCoef(ilpVariables_[*viBegin], *coefficient);             
    }
    ilpModel_.add(constraint);
}

template<class T>
template<class Iterator>
void
Cplex<T>::setStart(
    Iterator valueIterator
) {
    // delete existing mip starts
    if(ilpSolver_.getNMIPStarts() != 0) {
        ilpSolver_.deleteMIPStarts(0, ilpSolver_.getNMIPStarts());
    }

    // add new mip start
    for(size_t j=0; j<ilpVariables_.getSize(); ++j, ++valueIterator) {
        ilpStartValues_[j] = *valueIterator;
    }

    ilpSolver_.addMIPStart(ilpVariables_, ilpStartValues_);
}

template<class T>
inline void
Cplex<T>::setPreSolver(
    const PreSolver preSolver,
    const int passes 
) {
    switch(preSolver) {
    case PRE_SOLVER_NONE:
        ilpSolver_.setParam(IloCplex::PreInd, CPX_OFF);
        ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_OFF);
        return;
    case PRE_SOLVER_AUTO:
        ilpSolver_.setParam(IloCplex::PreInd, CPX_ON);
        ilpSolver_.setParam(IloCplex::PreDual, 0); // default
        ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_ON);
        ilpSolver_.setParam(IloCplex::RelaxPreDual, 0); // default
        break;
    case PRE_SOLVER_PRIMAL:
        ilpSolver_.setParam(IloCplex::PreInd, CPX_ON);
        ilpSolver_.setParam(IloCplex::PreDual, -1); 
        ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_ON);
        ilpSolver_.setParam(IloCplex::RelaxPreDual, -1);
        break;
    case PRE_SOLVER_DUAL:
        ilpSolver_.setParam(IloCplex::PreInd, CPX_ON);
        ilpSolver_.setParam(IloCplex::PreDual, 1); 
        ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_ON);
        ilpSolver_.setParam(IloCplex::RelaxPreDual, 1);
        break;
    }
}

template<class T>
inline void
Cplex<T>::setLPSolver(
    const LPSolver lpSolver
) {
    switch(lpSolver) {
    case LP_SOLVER_AUTO:
        ilpSolver_.setParam(IloCplex::RootAlg, 0); // default
        ilpSolver_.setParam(IloCplex::NodeAlg, 0); // default
        break;
    case LP_SOLVER_PRIMAL_SIMPLEX:
        ilpSolver_.setParam(IloCplex::RootAlg, 1);
        ilpSolver_.setParam(IloCplex::NodeAlg, 1);
        break;
    case LP_SOLVER_DUAL_SIMPLEX:
        ilpSolver_.setParam(IloCplex::RootAlg, 2);
        ilpSolver_.setParam(IloCplex::NodeAlg, 2);
        break;
    case LP_SOLVER_NETWORK_SIMPLEX:
        ilpSolver_.setParam(IloCplex::RootAlg, 3);
        ilpSolver_.setParam(IloCplex::NodeAlg, 3);
        break;
    case LP_SOLVER_BARRIER:
        ilpSolver_.setParam(IloCplex::RootAlg, 4);
        ilpSolver_.setParam(IloCplex::NodeAlg, 4);
        break;
    case LP_SOLVER_SIFTING:
        ilpSolver_.setParam(IloCplex::RootAlg, 5);
        ilpSolver_.setParam(IloCplex::NodeAlg, 5);
        break;
    case LP_SOLVER_CONCURRENT:
        ilpSolver_.setParam(IloCplex::RootAlg, 6);
        ilpSolver_.setParam(IloCplex::NodeAlg, 6);
        break;
    }
}

} // namespace ilp
} // namespace andres

#endif // #ifndef ANDRES_CPLEX_HXX
