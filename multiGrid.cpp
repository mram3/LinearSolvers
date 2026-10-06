/*
Compiling Instruction: g++ multiGrid.cpp src/*.cpp -Iinclude -std=c++17 -O3 -march=native -o main && ./main
*/

#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>

#include "LinearSolvers.h"

using namespace std;
using namespace std::chrono;

int main(){

    //Checking the multigrid method produces right answer by checking with direct methods
    int Nx = 48; //using coarse grid so that direct method can handle it
    int Ny = 48;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    int nCells = Nx*Ny;

    int nLevels = 3;
    int divisor = pow(2, nLevels-1);

    if(Nx % divisor != 0 || Ny % divisor != 0){
        cout << "Cell counts must be divisible by " << divisor <<" to have "<< nLevels<< " levels\n";
        return 1;
    }
    
    vector<Grid> mesh; //creating a vector of grids to store multi level meshes
    mesh.reserve(nLevels);
    vector<Equations> levels; //creating a vector of equations to store matrices at multi lvl
    levels.reserve(nLevels);

    for(int i = 0; i < nLevels; ++i){
        mesh.emplace_back(Nx, Ny, xmin, xmax, ymin, ymax);

        levels.emplace_back(&mesh[i]);

        levels[i].assemblePoissonMatrix();

        Nx /= 2; Ny /= 2;
    }

    LinearSolvers solve;

    //Getting the exact solution using exact LU
    std::vector<double> x_exact(nCells, 0.0);
    
    auto start = high_resolution_clock::now();
    solve.directLU(levels[0].A, levels[0].b, x_exact);
    auto end = high_resolution_clock::now();
    double timeLU = duration<double, std::milli>(end - start).count();

    vector<double> x(nCells, 0.0); //Solution through multigrid methods

    start = high_resolution_clock::now();
    int iterV = solve.GMG(levels, x, 1e-6, Cycle::vCycle);
    end = high_resolution_clock::now();
    double timeV = duration<double, milli>(end - start).count();
    double errV = MathTools::L2Norm(MathTools::vectorSub(x, x_exact));

    fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterW = solve.GMG(levels, x, 1e-6, Cycle::wCycle);
    end = high_resolution_clock::now();
    double timeW = duration<double, milli>(end - start).count();
    double errW = MathTools::L2Norm(MathTools::vectorSub(x, x_exact));

    fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterF = solve.GMG(levels, x, 1e-6, Cycle::fCycle);
    end = high_resolution_clock::now();
    double timeF = duration<double, milli>(end - start).count();
    double errF = MathTools::L2Norm(MathTools::vectorSub(x, x_exact));

    //printing results to confirm solution is same

    cout << "\n===============================================================================\n";
    cout << left << setw(12) << "Cycle" 
                 << setw(17) << "No. of Cycles"
                 << setw(17) << "Wall Time (ms)"
                 << "Error Vs Exact LU" << endl;
    cout << "================================================================================\n";

    cout << scientific << setprecision(4);

    cout << left << setw(12) << "V Cycle" << setw(17) << iterV << fixed << setw(17) << timeV << scientific << errV << endl;
    cout << left << setw(12) << "W Cycle" << setw(17) << iterW << fixed << setw(17) << timeW << scientific << errW << endl;
    cout << left << setw(12) << "F Cycle" << setw(17) << iterF << fixed << setw(17) << timeF << scientific << errF << endl;
    cout << "--------------------------------------------------------------------------------\n";

    cout << left << setw(12) << "Exact LU" << setw(17) << "N/A" << fixed << setw(17) << timeLU << scientific << "0.000" << endl;
    cout << "================================================================================\n";
    //Now moving to much more finer mesh
    Nx = 2048;
    Ny = 2048;

    nCells = Nx*Ny;

    nLevels = 8;
    divisor = pow(2, nLevels-1);

    if(Nx % divisor != 0 || Ny % divisor != 0){
        cout << "Cell counts must be divisible by " << divisor <<" to have "<< nLevels<< " levels\n";
        return 1;
    }

    mesh.clear();
    levels.clear();
    x.clear();

    mesh.reserve(nLevels);
    levels.reserve(nLevels);

    for(int i = 0; i < nLevels; ++i){
        mesh.emplace_back(Nx, Ny, xmin, xmax, ymin, ymax);

        levels.emplace_back(&mesh[i]);

        levels[i].assemblePoissonMatrix();

        Nx /= 2; Ny /= 2;
    }
    
    x.assign(nCells, 0.0);
    start = high_resolution_clock::now();
    iterV = solve.GMG(levels, x, 1e-6, Cycle::vCycle);
    end = high_resolution_clock::now();
    timeV = duration<double, milli>(end - start).count();

    x.assign(nCells, 0.0);
    start = high_resolution_clock::now();
    iterW = solve.GMG(levels, x, 1e-6, Cycle::wCycle);
    end = high_resolution_clock::now();
    timeW = duration<double, milli>(end - start).count();

    x.assign(nCells, 0.0);
    start = high_resolution_clock::now();
    iterF = solve.GMG(levels, x, 1e-6, Cycle::fCycle);
    end = high_resolution_clock::now();
    timeF = duration<double, milli>(end - start).count();


    cout << "\n===============================================================================\n";
    cout << left << setw(12) << "Cycle" 
                 << setw(17) << "No. of Cycles"
                 << setw(17) << "Wall Time (ms)" << endl;
    cout << "================================================================================\n";

    cout << left << setw(12) << "V Cycle" << setw(17) << iterV << fixed << setw(17) << timeV << endl;
    cout << left << setw(12) << "W Cycle" << setw(17) << iterW << fixed << setw(17) << timeW << endl;
    cout << left << setw(12) << "F Cycle" << setw(17) << iterF << fixed << setw(17) << timeF << endl;
    cout << "--------------------------------------------------------------------------------\n";

    return 0;
}