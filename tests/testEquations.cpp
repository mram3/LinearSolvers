/*
Compiling Instructions: g++ tests/testEquations.cpp src/*.cpp -Iinclude -std=c++17 -o test && ./test
*/

#include "Equations.h"

#include <iostream>

using namespace std;

int main(){
    int nx = 3, ny = 2, nCells = nx*ny;
    Grid mesh(nx, ny, 0.0, 1.5, 0.0, 1.0);
    Equations eqn;

    eqn.assemblePoissonMatrix(mesh);

    //Unpacking CSR into a Dense Matrix
    std::vector<std::vector<double>> M(nx*ny, std::vector<double>(nx*ny, 0.0));
    const auto& rowPtr = eqn.A.getrowPtr();
    const auto& col = eqn.A.getcol();
    const auto& values = eqn.A.getvalues();

    for(int i = 0; i < rowPtr.size()-1; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            M[i][col[j]] = values[j];
        }
    }

    for(int i = 0; i < nCells; ++i){
        for(int j = 0; j< nCells; ++j){
            cout << M[i][j] << "  ";
        }
        cout << endl;
    }

    cout << endl;

    for(int i = 0; i < nCells; ++i){
        cout << eqn.b[i] << endl;
    }

    return 0;
}