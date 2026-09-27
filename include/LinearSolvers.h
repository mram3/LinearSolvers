#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "Preconditioners.h"
#include "Equations.h"
#include "MathTools.h"
#include "Type.h"

class LinearSolvers{
public:

    int richardsonIteration
    (
        const Matrix& A,
        const Preconditioners& M,
        const std::vector<double>& b,
        std::vector<double>& x,
        double tolerance,
        Side side
    );

    void directCholesky
    (
        const Matrix& A_sparse, 
        const std::vector<double>& b, 
        std::vector<double>& x
    );

    static void directLU
    (
        const Matrix& A_sparse, 
        const std::vector<double>& b, 
        std::vector<double>& x
    );

    int GMG
    (
        std::vector<Equations>& levels,
        std::vector<double>& x,
        double tolarance,
        Cycle cycle
    );

    static void gaussSeidel
    (
        const Matrix& A, 
        const std::vector<double>& b,
        std::vector<double>& x,
        int maxIter
    );

};
#endif
