#include "Grid.h"

Grid::Grid
(
    int Nx_, int Ny_,
    double xmin_, double xmax_, double ymin_, double ymax_
) : Nx(Nx_), Ny(Ny_), xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_)
{
    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;

    x.reserve(Nx);
    y.reserve(Ny);
    
    for(int i = 0; i < Nx; ++i){
        x[i] = xmin + (i+0.5) * dx;
    }

    for(int j = 0; j < Ny; ++j){
        y[j] = ymin + (j+0.5) * dy;
    }
}