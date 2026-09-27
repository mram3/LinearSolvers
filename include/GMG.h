#pragma once

#include "Equations.h"

class GMG
{
public:
    static void vCycle
    (
        std::vector<Equations>& level,
        std::vector<double>& result
    );

    static void wCycle
    (
        std::vector<Equations>& level,
        std::vector<double>& result
    );
    static void fCycle
    (
        std::vector<Equations>& level,
        std::vector<double>& result
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