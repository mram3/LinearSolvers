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
    int xCells =  mesh->Nx, yCells = mesh->Ny, rowIdx;

    double dx = mesh->dx, dy = mesh->dy;
    double xmin = mesh->xmin, xmax = mesh->xmax, ymin = mesh->ymin, ymax = mesh->ymax;

    //Physical Properties 

    double gamma = 0.5; //diffusion coefficient
    double rho = 1.0; //density
    double u = 5.0;   //x velocity
    double v = 5.0;   //y velocity

    //Boundary Conditions

    double phi_east = 0.0;
    double phi_west = 0.0;
    double phi_north = 0.0;
    double phi_south = 0.0;
    double phi_south_ = 100.0;

    //mass flow rates
    double mW = rho*u*dy, mE = rho*u*dy, mN = rho*v*dx, mS = rho*v*dx;

    //rhs
    b.assign(xCells * yCells, 0.0);

    for(int j = 0; j < yCells; ++j){
        for(int i = 0; i < xCells; ++i){
            //coefficients
            double aW = -gamma*dy/dx,aE = -gamma*dy/dx,aN = -gamma*dx/dy,aS = -gamma*dx/dy;
            double aP = 0.0, Su = 0.0, Sp =0.0;
            
            if(i == 0){
                aW = 0.0;
                aE -= max(-mE, 0.0);
                Su += phi_west * max(mW, 0.0);
                Su += 2*gamma*dy*phi_west/dx;
                Sp += 2*gamma*dy/dx;
                Sp += max(mW, 0.0);
            }
            else if (i == xCells-1){
                aE = 0.0;
                aW -= max(mW, 0.0);
                Su += phi_east * max(-mE, 0.0);
                Su += 2*gamma*dy*phi_east/dx;
                Sp += 2*gamma*dy/dx;
                Sp += max(-mE, 0.0);
            }
            else{
                aW -= max(mW, 0.0);
                aE -= max(-mE, 0.0);
            }

            if(j == 0){
                aS = 0.0;
                aN -= max(-mN, 0.0);
                Su += (mesh->x[i] > 0.4 && mesh->x[i] < 0.6) 
                                                ? phi_south_*(2*gamma*dx/dy+max(mS, 0.0))
                                                : phi_south*(2*gamma*dx/dy+max(mS, 0.0));
                Sp += 2*gamma*dx/dy;
                Sp += max(mS, 0.0);
            }
            else if(j == yCells - 1){
                aN = 0.0;
                aS -= max(mS, 0.0);
                Su += phi_north * max(-mN, 0.0);
                Su += 2*gamma*dx*phi_north/dy;
                Sp += 2*gamma*dx/dy;
                Sp += max(-mN, 0.0);
            }
            else{
                aN -= max(-mN, 0.0);
                aS -= max(mS, 0.0);
            }
            
            aP = -(aW + aE + aN + aS - Sp);

            rowIdx = j * xCells + i;

            b[rowIdx] = Su;

            A.addCoeff(rowIdx, rowIdx, aP);

            if(aW != 0.0) A.addCoeff(rowIdx, rowIdx - 1, aW);            
            if(aE != 0.0) A.addCoeff(rowIdx, rowIdx + 1, aE);            
            if(aN != 0.0) A.addCoeff(rowIdx, rowIdx + xCells, aN);            
            if(aS != 0.0) A.addCoeff(rowIdx, rowIdx - xCells, aS);            
        }
    }
    
    A.compressToCSR();
}