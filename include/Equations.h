#ifndef EQUATIONS_H
#define EQUATIONS_H

#include "Grid.h"
#include "Matrix.h"

#include <vector>

class Equations{
public:

    const Grid* mesh = nullptr;
    Matrix A;
    std::vector<double> b;

    //default constructor
    Equations() = default;

    //Parameterized Constructor
    Equations(const Grid* mesh_) : mesh(mesh_) {};

    void assemblePoissonMatrix();

    void assembleConvectionDiffusionEquation();

};

#endif