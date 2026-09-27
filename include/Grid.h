#ifndef GRID_H
#define GRID_H

#include <vector>

class Grid 
{
public:
    int Nx, Ny;
    double xmin, xmax, ymin, ymax, dx, dy;
    std::vector<double> x, y;

    Grid
    (
        int Nx, int Ny, 
        double xmin, double xmax, double ymin, double ymax
    );

};

#endif