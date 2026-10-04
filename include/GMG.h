#pragma once

#include "Equations.h"
#include "Type.h"
#include "Smoothers.h"
#include "MathTools.h"

class GMG
{
public:
    static void runCycles
    (
        int lvl,
        std::vector<Equations>& levels,
        std::vector<double>& x,
        const std::vector<double>& rhs,
        Cycle cycle
    );

private:
    static std::vector<double> restrictRes
    (
        const Grid* meshFine,
        const Grid* meshCoarse,
        const std::vector<double>& prevRes
    );

    static std::vector<double> prolongErr
    (
        const Grid* meshCoarse,
        const Grid* meshFine,
        const std::vector<double>& prevErr
    );
};