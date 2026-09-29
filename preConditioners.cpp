/*
Compiling Instruction: g++ preConditioners.cpp src/*.cpp -Iinclude -std=c++17 -O3 -march=native -o main && ./main
*/
#include <iostream>
#include <vector>
#include <algorithm> 
#include <chrono>  
#include <iomanip>   

#include "LinearSolvers.h"

using namespace std;
using namespace std::chrono;

int main(){
    int Nx = 10;
    int Ny = 10;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    int NCells = Nx*Ny;

    std::vector<double> x(NCells, 0.0); // initialize the solution vector

    Grid mesh(Nx, Ny, xmin, xmax, ymin, ymax);
    Equations eqn(&mesh); 
    eqn.assemblePoissonMatrix();
    LinearSolvers solve;
 
    //eqn.assemblePoissonMatrix(Nx, Ny); // creates a square matrix of NCells*NCells
    
    jacobiPreconditioner jP; 
    sorPreconditioner gS(1.0);
    sorPreconditioner sP(1.56);
    sgsPreconditioner sGS;
    ilu0Preconditioner ilu;
    ic0Preconditioner ic;

    //setting up the preconditioners
    jP.setup(eqn.A);
    gS.setup(eqn.A);
    sP.setup(eqn.A);
    sGS.setup(eqn.A);
    ilu.setup(eqn.A);
    ic.setup(eqn.A);

    double timeJP, timeGS, timeSP, timeSGS, timeILU, timeIC, timeLU, timeChol;
    double errJP, errGS, errSP, errSGS, errILU, errIC, errChol;

    //Getting the exact solution using exact LU
    std::vector<double> x_exact(NCells, 0.0);
    
    auto start = high_resolution_clock::now();
    solve.directLU(eqn.A, eqn.b, x_exact);
    auto end = high_resolution_clock::now();
    timeLU = duration<double, std::milli>(end - start).count();

    //Iterating the preconditioned system

    // 1. Jacobi Preconditioner
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterJP = solve.richardsonIteration(eqn.A, jP, eqn.b, x, 1e-6, Side::Left);
    end = high_resolution_clock::now();
    timeJP = duration<double, std::milli>(end - start).count();
    errJP = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    // 2. Gauss-Seidel
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterGS = solve.richardsonIteration(eqn.A, gS, eqn.b, x, 1e-6, Side::Left);
    end = high_resolution_clock::now();
    timeGS = duration<double, std::milli>(end - start).count();
    errGS = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    // 3. SOR
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterSP = solve.richardsonIteration(eqn.A, sP, eqn.b, x, 1e-6, Side::Left);
    end = high_resolution_clock::now();
    timeSP = duration<double, std::milli>(end - start).count();
    errSP = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    // 4. Symmetric GS
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterSGS = solve.richardsonIteration(eqn.A, sGS, eqn.b, x, 1e-6, Side::Left);
    end = high_resolution_clock::now();
    timeSGS = duration<double, std::milli>(end - start).count();
    errSGS = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    // 5. ILU(0)
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterILU = solve.richardsonIteration(eqn.A, ilu, eqn.b, x, 1e-6, Side::Left);
    end = high_resolution_clock::now();
    timeILU = duration<double, std::milli>(end - start).count();
    errILU = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    // 6. IC(0)
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    int iterIC = solve.richardsonIteration(eqn.A, ic, eqn.b, x, 1e-6, Side::Left);
    end = high_resolution_clock::now();
    timeIC = duration<double, std::milli>(end - start).count();
    errIC = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    // 7. Exact Cholesky
    std::fill(x.begin(), x.end(), 0.0);
    start = high_resolution_clock::now();
    solve.directCholesky(eqn.A, eqn.b, x);
    end = high_resolution_clock::now();
    timeChol = duration<double, std::milli>(end - start).count();
    errChol = MathTools::L2Norm(MathTools::vectorSub(x_exact, x));

    //Printing results
    
    cout << "\n======================================================================================\n";
    cout << left << setw(30) << "Solver / Preconditioner" 
         << setw(15) << "Iterations" 
         << setw(20) << "Wall Time (ms)" 
         << "Error vs Exact LU" << endl;
    cout << "======================================================================================\n";
    
    cout << scientific << setprecision(4);
    
    cout << left << setw(30) << "Jacobi (Left)"      << setw(15) << iterJP  << fixed << setprecision(3) << setw(20) << timeJP << scientific << errJP << endl;
    cout << left << setw(30) << "Gauss-Seidel (Left)"<< setw(15) << iterGS  << fixed << setw(20) << timeGS << scientific << errGS << endl;
    cout << left << setw(30) << "SOR w=1.56 (Left)"  << setw(15) << iterSP  << fixed << setw(20) << timeSP << scientific << errSP << endl;
    cout << left << setw(30) << "Symmetric GS (Left)"<< setw(15) << iterSGS << fixed << setw(20) << timeSGS << scientific << errSGS << endl;
    cout << left << setw(30) << "ILU(0) (Left)"      << setw(15) << iterILU << fixed << setw(20) << timeILU << scientific << errILU << endl;
    cout << left << setw(30) << "IC(0) (Left)"       << setw(15) << iterIC  << fixed << setw(20) << timeIC << scientific << errIC << endl;
    
    cout << "--------------------------------------------------------------------------------------\n";
    
    //Error for Exact LU is 0.0 because it is the reference vector itself
    cout << left << setw(30) << "Exact Dense LU"       << setw(15) << "N/A" << fixed << setw(20) << timeLU   << scientific << 0.0     << endl;
    cout << left << setw(30) << "Exact Dense Cholesky" << setw(15) << "N/A" << fixed << setw(20) << timeChol << scientific << errChol << endl;
    cout << "======================================================================================\n";

    return 0;
}