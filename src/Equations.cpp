#include "Equations.h"

#include <cmath>

using namespace std;

const double PI = acos(-1.0);

void Equations::assemblePoissonMatrix() 
{
    int xCells = mesh->Nx, yCells = mesh->Ny, rowIdx;
    double dx = mesh->dx, dy = mesh->dy, dV = dx*dy;

    //Boundary conditions
    double phiEast = 0.0;
    double phiWest = 0.0;
    double phiNorth = 0.0;
    double phiSouth = 0.0;
    
    b.assign(xCells * yCells, 0.0);

    for(int j = 0; j < yCells; ++j){
        for(int i = 0; i < xCells; ++i){

            double aW = -1.0*dy/dx, aE = -1.0*dy/dx, aS = -1.0*dx/dy, aN = -1.0*dx/dy, aP;
            double Su = 0.0, Sp = 0.0;

            if(i == 0){
                aW = 0.0;
                Su += phiWest * 2.0*dy/dx;
                Sp += 2.0*dy/dx;
            }
            else if(i == xCells-1){
                aE = 0.0;
                Su += phiEast * 2.0*dy/dx;
                Sp += 2.0*dy/dx;
            }

            if(j == 0){
                aS = 0.0;
                Su += phiSouth * 2.0*dx/dy;
                Sp += 2.0*dy/dx;
            }
            else if(j == yCells-1){
                aN = 0.0;
                Su += phiNorth * 2.0*dx/dy;
                Sp += 2.0*dx/dy;
            }

            aP = - (aW+aE+aN+aS-Sp);

            rowIdx = j * xCells + i;
            A.addCoeff(rowIdx, rowIdx, aP);

            if(aW != 0.0) A.addCoeff(rowIdx, rowIdx-1, aW);
            if(aE != 0.0) A.addCoeff(rowIdx, rowIdx+1, aE);
            if(aS != 0.0) A.addCoeff(rowIdx, rowIdx-xCells,aS);
            if(aN != 0.0) A.addCoeff(rowIdx, rowIdx+xCells,aN);


            b[rowIdx] = 2*PI*PI*sin(PI*mesh->x[i])*sin(PI*mesh->y[j])*dV + Su;
        }
    }

    //Comprssing COO to CSR
    A.compressToCSR();
}

void Equations::assembleConvectionDiffusionEquation()
{
    int Nx =  mesh->Nx, Ny = mesh->Ny, rowIdx;

    double dx = mesh->dx, dy = mesh->dy;
    double xmin = mesh->xmin, xmax = mesh->xmax, ymin = mesh->ymin, ymax = mesh->ymax;

    
}