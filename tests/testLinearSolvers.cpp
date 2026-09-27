/*
Compiling Instructions: g++ tests/testLinearSolvers.cpp src/*.cpp -Iinclude -std=c++17 -o test && ./test
*/

#include "LinearSolvers.h"

#include <iostream>

using namespace std;

int main(){
    int nx = 3, ny = 2, nCells = nx*ny;
    Grid mesh(nx, ny, 0.0, 1.5, 0.0, 1.0);
    Equations eqn;

    eqn.assemblePoissonMatrix(mesh);

    vector<double> x(nCells, 0.0), xLU(nCells, 0.0);
    LinearSolvers::gaussSeidel(eqn.A, eqn.b, x, 2000);
    LinearSolvers::directLU(eqn.A, eqn.b, xLU);

    for(int i = 0; i < nCells; ++i){
        cout << x[i] << " " << xLU[i] << endl;
    }

    return 0;
}