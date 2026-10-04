#include "GMG.h"

using namespace std;

std::vector<double> GMG::restrictRes
(
    const Grid* meshFine,
    const Grid* meshCoarse,
    const vector<double>& prevRes
)
{
    int Nx = meshCoarse->Nx, Ny = meshCoarse->Ny, currResIdx, prevResIdx;
    int nx = meshFine->Nx;

    vector<double> currRes(Nx*Ny, 0.0);

    for(int jC = 0;  jC < Ny; ++jC){
        for(int iC = 0; iC < Nx; ++iC){
            currResIdx = jC * Nx + iC;

            int iF = 2*iC;//coords in fine mesh
            int jF = 2*jC;

            prevResIdx = jF * nx + iF;//bottom left cell

            currRes[currResIdx]= prevRes[prevResIdx] //bottom left cell
                                +prevRes[prevResIdx+1] //bottom right cell
                                +prevRes[prevResIdx+nx] //top left cell
                                +prevRes[prevResIdx+nx+1]; //top right cell
        }
    }
    return currRes;
}

std::vector<double> GMG::prolongErr
(
    const Grid* meshCoarse,
    const Grid* meshFine,
    const std::vector<double>& prevErr
)
{
    int Nx = meshCoarse->Nx, Ny = meshCoarse->Ny, currErrIdx, prevErrIdx;
    int nx = meshFine->Nx, ny = meshFine->Ny;
    vector<double> currErr(nx*ny, 0.0);

    for(int jC = 0; jC < Ny; ++jC){
        for(int iC = 0; iC < Nx; ++iC){
            prevErrIdx = jC * Nx + iC;

            int iF = 2*iC; //coords in fine mesh
            int jF = 2*jC;

            currErrIdx = jF * nx + iF; //bottom left fine mesh

            //zero order prolongation
            //error on the parent Coarse cell is passed to fine cells
            currErr[currErrIdx] = prevErr[prevErrIdx];
            currErr[currErrIdx+1] = prevErr[prevErrIdx]; //bottom right
            currErr[currErrIdx + nx] = prevErr[prevErrIdx]; //top left
            currErr[currErrIdx + nx +1] = prevErr[prevErrIdx]; //top right
        }
    }
    return currErr;
}

void GMG::runCycles
(
    int lvl,
    vector<Equations>& levels,
    vector<double>& x,
    const vector<double>& rhs,
    Cycle cycle
)
{
    int n = levels.size();
    if(lvl == n-1){
        Smoothers::gaussSeidel(levels[lvl].A, rhs, x, 50);
        return;
    }

    //Pre Sweeps
    //Few Gauss iterations on finer grids
    Smoothers::gaussSeidel(levels[lvl].A, rhs, x, 3);

    //calculating the residual at current level
    vector<double> currRes = MathTools::vectorSub(rhs, levels[lvl].A.SpMV(x));

    //restricting the residual to next coarse level
    vector<double> nextRes = restrictRes(levels[lvl].mesh, levels[lvl+1].mesh, currRes);

    //initiating error array for next coarse level
    vector<double> error(levels[lvl+1].b.size());

    switch (cycle){
        case Cycle::vCycle:
            runCycles(lvl+1, levels, error, nextRes, Cycle::vCycle);
            break;
        
        case Cycle::wCycle:
            runCycles(lvl+1, levels, error, nextRes, Cycle::wCycle);
            runCycles(lvl+1, levels, error, nextRes, Cycle::wCycle);
            break;
        
        case Cycle::fCycle:
            runCycles(lvl+1, levels, error, nextRes, Cycle::fCycle);
            runCycles(lvl+1, levels, error, nextRes, Cycle::vCycle);
            break;
    }

    vector<double> errCorr = prolongErr(levels[lvl+1].mesh, levels[lvl].mesh, error);
    x = MathTools::vectorAdd(x, errCorr);
    Smoothers::gaussSeidel(levels[lvl].A, rhs, x, 3);
}

/*void GMG::vCycle
(
    vector<Equations>& levels,
    vector<double>& x
)
{
    int n = levels[0].b.size(), nLevels = levels.size();
    vector<vector<double>> residuals(nLevels);
    vector<vector<double>> error(nLevels-1);

    //1. Pre Sweeps
    //Few Gauss Seidel iteration on fine grid to eliminate high frequency errors
    Smoothers::gaussSeidel(levels[0].A, levels[0].b, x, 3); //smoother

    //2. Restirction (Going down to coarse grids)
    residuals[0] = MathTools::vectorSub(levels[0].b, levels[0].A.SpMV(x));

    for(int i = 0; i < nLevels-1; ++i){

        residuals[i+1] = restrictRes(levels[i].mesh, levels[i+1].mesh, residuals[i]);

        error[i].assign(levels[i+1].b.size(), 0.0);

        if(i == nLevels - 2){
            Smoothers::gaussSeidel(levels[i+1].A, residuals[i+1], error[i], 50);
            continue;
        }

        Smoothers::gaussSeidel(levels[i+1].A, residuals[i+1], error[i], 3);

        residuals[i+1] = MathTools::vectorSub(residuals[i+1], levels[i+1].A.SpMV(error[i]));
    }

    //3. prolongation (Going up to fine grids)
    for(int j = nLevels-2; j >= 0; --j){

        auto errCorr = prolongErr(levels[j+1].mesh, levels[j].mesh, error[j]);

        if(j == 0){
            x = MathTools::vectorAdd(x, errCorr);
            break;
        }
        
        error[j-1] = MathTools::vectorAdd(error[j-1], errCorr);

        Smoothers::gaussSeidel(levels[j].A, residuals[j], error[j-1], 3);
    }

    //4.Correction and final iterations 
    //Post Sweeps
    Smoothers::gaussSeidel(levels[0].A, levels[0].b, x, 3);
}*/