#pragma once

#include "Equations.h"

class Smoothers{
public:
    static void gaussSeidel
    (
        Matrix& A,
        const std::vector<double>& rhs,
        std::vector<double>& x, 
        int maxIter
    );
};