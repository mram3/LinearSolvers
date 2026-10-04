/*
Compiling Instruction: g++ multiGrid.cpp src/*.cpp -Iinclude -std=c++17 -O3 -march=native -o main && ./main
*/

#include <iostream>
#include <cmath>
#include <chrono>

#include "LinearSolvers.h"

using namespace std;
using namespace std::chrono;

int main(){
    int Nx = 4096;
    int Ny = 4096;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    int nCells = Nx*Ny;

    int nLevels = 8;
    int divisor = pow(2, nLevels-1);

    if(Nx % divisor != 0 || Ny % divisor != 0){
        cout << "Cell counts must be divisible by " << divisor <<" to have "<< nLevels<< " levels\n";
        return 1;
    }
    
    vector<Grid> mesh;
    mesh.reserve(nLevels);
    vector<Equations> levels;
    levels.reserve(nLevels);

    for(int i = 0; i < nLevels; ++i){
        mesh.emplace_back(Nx, Ny, xmin, xmax, ymin, ymax);

        levels.emplace_back(&mesh[i]);

        levels[i].assemblePoissonMatrix();

        Nx /= 2; Ny /= 2;
    }

    vector<double> x(nCells, 0.0);
    LinearSolvers solve;

    auto start = high_resolution_clock::now();
    int iterV = solve.GMG(levels, x, 1e-6, Cycle::vCycle);
    auto end = high_resolution_clock::now();
    double timeV = duration<double, milli>(end - start).count();

    fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterW = solve.GMG(levels, x, 1e-6, Cycle::wCycle);
    end = high_resolution_clock::now();
    double timeW = duration<double, milli>(end - start).count();

    fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterF = solve.GMG(levels, x, 1e-6, Cycle::fCycle);
    end = high_resolution_clock::now();
    double timeF = duration<double, milli>(end - start).count();

    cout << iterV << " " << timeV << endl;
    cout << iterW << " " << timeW << endl;
    cout << iterF << " " << timeF << endl;
}