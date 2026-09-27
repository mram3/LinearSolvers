#include "LinearSolvers.h"
#include "GMG.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

int LinearSolvers::richardsonIteration
(
    const Matrix& A,
    const Preconditioners& M,
    const vector<double>& b,
    vector<double>& x,
    double tolrance,
    Side side
)
{
    int max_iter = 10000;
    int iter = 0;
    double residual = 1.0;
    vector<double> z(b.size(), 0.0);//initializing the M^-1 * (b - Ax) vector 
    vector<double> z_lower(b.size(), 0.0);

    while(iter < max_iter){

        //A*x
        vector<double> Ax = A.SpMV(x);
        //b-A*x
        vector<double> r = MathTools::vectorSub(b, Ax);

        if(residual < tolrance){
            return iter;
        }
        
        switch (side)
        {
        case Side::Left:
            M.apply(z, r);
            residual = MathTools::L2Norm(r);
            x = MathTools::vectorAdd(z, x);
            break;

        case Side::Right:
            residual = MathTools::L2Norm(r);
            M.apply(z, r);
            x = MathTools::vectorAdd(z, x);
            break;

        case Side::Split:
            
            // 1. Forward sweep only: L * z_lower = r
            M.applyLower(z_lower, r); 
            
            residual = MathTools::L2Norm(z_lower);
            
            // 2. Backward sweep only: U * z = z_lower
            M.applyUpper(z, z_lower); 
            
            x = MathTools::vectorAdd(z, x);
            break;

        case Side::None:
            residual = MathTools::L2Norm(r);
            x = MathTools::vectorAdd(r, x);
            break;
        }

        if(residual < tolrance){
            return iter;
        }
        
        iter++;
    }

    cout << "Does not converge\n";
    return iter;
}

void LinearSolvers::directCholesky
(
    const Matrix& A_sparse,
    const std::vector<double>& b,
    std::vector<double>& x
)
{
    int n = b.size();
    
    //Unpacking CSR into a Dense Matrix
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    const auto& rowPtr = A_sparse.getrowPtr();
    const auto& col = A_sparse.getcol();
    const auto& values = A_sparse.getvalues();

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            A[i][col[j]] = values[j];
        }
    }

    //Exact Dense Cholesky Factorization (A = L * L^T)
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    
    for(int i = 0; i < n; ++i){
        for(int j = 0; j <= i; ++j){
            double sum = 0.0;
            
            // Calculating the dot product of previous elements
            for(int k = 0; k < j; ++k){
                sum += L[i][k] * L[j][k];
            }
            
            // Applying exact Cholesky equations
            if(i == j){
                L[i][j] = std::sqrt(A[i][i] - sum);
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }

    //Forward Substitution (L * y = b)
    std::vector<double> y(n, 0.0);
    for(int i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = 0; j < i; ++j){
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    //Backward Substitution (L^T * x = y)
    for(int i = n - 1; i >= 0; --i){
        double sum = 0.0;
        for(int j = i + 1; j < n; ++j){
            sum += L[j][i] * x[j]; 
        }
        x[i] = (y[i] - sum) / L[i][i];
    }
}

void LinearSolvers::directLU
(
    const Matrix& A_sparse,
    const std::vector<double>& b,
    std::vector<double>& x
)
{
    int n = b.size();

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    const auto& rowPtr = A_sparse.getrowPtr();
    const auto& col = A_sparse.getcol();
    const auto& values = A_sparse.getvalues();

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            A[i][col[j]] = values[j];
        }
    }

    for(int i = 0; i < n; ++i){
        
        for(int j = i; j < n; ++j){
            double sum = 0.0;
            for(int k = 0; k < i; ++k){
                sum += A[i][k] * A[k][j]; // L_ik * U_kj
            }
            A[i][j] = A[i][j] - sum;
        }

        for(int j = i + 1; j < n; ++j){
            double sum = 0.0;
            for(int k = 0; k < i; ++k){
                sum += A[j][k] * A[k][i]; // L_jk * U_ki
            }
            A[j][i] = (A[j][i] - sum) / A[i][i]; // Divide by U_ii
        }
    }

    //Forward Substitution (L * y = b)
    std::vector<double> y(n, 0.0);
    for(int i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = 0; j < i; ++j){
            sum += A[i][j] * y[j]; // A[i][j] acts as L
        }
        y[i] = b[i] - sum;
    }

    //Backward Substitution (U * x = y)
    for(int i = n - 1; i >= 0; --i){
        double sum = 0.0;
        for(int j = i + 1; j < n; ++j){
            sum += A[i][j] * x[j]; // A[i][j] acts as U
        }
        x[i] = (y[i] - sum) / A[i][i]; // A[i][i] acts as U_ii
    }
}

int LinearSolvers::GMG
(
    std::vector<Equations>& levels,
    std::vector<double>& x,
    double tolarance,
    Cycle cycle
)
{
    double r, r0;
    int iter = 0, maxIter = 1000;

    vector<double> z = levels[0].A.SpMV(x);
    r0 = MathTools::L2Norm(MathTools::vectorSub(levels[0].b, z));
    r = r0;

    if(r0 < std::numeric_limits<double>::epsilon()) return iter;

    switch (cycle)
    {
    case Cycle::vCycle:
        while(iter < maxIter && r/r0 > tolarance){
            GMG::vCycle(levels, x);
            z = levels[0].A.SpMV(x);
            r = MathTools::L2Norm(MathTools::vectorSub(levels[0].b, z));
            iter++;
        }
        if(iter==maxIter) cout << "Does not converge.\n";
        break;
    }

    return iter;
}

void LinearSolvers::gaussSeidel
(
    const Matrix& A, 
    const std::vector<double>& b,
    std::vector<double>& x,
    int maxIter
)
{
    const auto& rowPtr = A.getrowPtr();
    const auto& col = A.getcol();
    const auto& values = A.getvalues();

    int n = rowPtr.size() - 1, iter = 0;
    vector<double> diag(n, 1.0);

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                diag[i] = values[j];
                break;
            }
        }
    }

    while(iter < maxIter){
        for(int i = 0; i < n; ++i){
            double sum = b[i];
            for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
                if(i != col[j]){
                    sum -= values[j]*x[col[j]];
                }
            }
            x[i] = sum/diag[i];
        }
        iter++;
    }
}