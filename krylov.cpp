/*
Compiling Instruction: g++ krylov.cpp src/*.cpp -Iinclude -std=c++17 -O3 -march=native -o main && ./main
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>

#include "LinearSolvers.h"

using namespace std::chrono;

int main(){

    int Nx = 2048, Ny = 2048;
    int nCells = Nx * Ny;

    double xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0;
    Grid mesh(Nx, Ny, xmin, xmax, ymin, ymax);

    Equations eqn(&mesh);
    eqn.assembleConvectionDiffusionEquation();

    std::vector<double> x(nCells, 0.0);

    sgsPreconditioner M;
    M.setup(eqn.A);

    LinearSolvers solve;

    auto start = high_resolution_clock::now();
    int iter = solve.PBiCGStab(eqn.A, M, x, eqn.b, 1e-6);
    auto end = high_resolution_clock::now();
    double timeBCG = duration<double, std::milli>(end - start).count();

    std::cout << iter << "\t" << timeBCG << std::endl;

    return 0;
}
