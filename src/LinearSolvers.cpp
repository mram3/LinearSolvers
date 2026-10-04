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

    while(iter < maxIter && r/r0 > tolarance){
        
        //triggering the recursion here
        GMG::runCycles(0, levels, x, levels[0].b, cycle);

        //calculating the residual
        z = levels[0].A.SpMV(x);
        r = MathTools::L2Norm(MathTools::vectorSub(levels[0].b, z));
        iter++;
    }

    if(iter==maxIter) cout << "Does not converge.\n";
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

int LinearSolvers::PCG
(
    const Matrix& A,
    const Preconditioners& M,
    vector<double>& x,
    const vector<double>& rhs,
    double tolarance
)
{
    int maxIter = 10000, iter = 0;

    //brief PCG (Yousef Saad - Iterative Methods for Sparse Linear Systems)
    //computing r0 := b - A*x0
    vector<double> res = MathTools::vectorSub(rhs, A.SpMV(x));
    //z0 := inv(M) * r0
    vector<double> z(res.size());
    M.apply(z, res);

    //p0 := z0
    vector<double> searchDirn = z;

    double resz = MathTools::innerProd(z, res);
    double r0 = MathTools::L2Norm(res);
    double r = r0;
    
    if(r0 < std::numeric_limits<double>::epsilon()) return iter;
    
    while(iter < maxIter && r/r0 > tolarance){
        vector<double> Ap = A.SpMV(searchDirn);

        //alpha(j) := tr(r(j))*z(j) / (tr(p(j) * A * p(j))
        double alpha = resz/MathTools::innerProd(searchDirn, Ap);

        //x(j+1) := x(j) + alpha(j) * p(j)
        MathTools::daxpy(x, alpha, searchDirn);

        //r(j+1) := r(j) - alpha(j)*A*p(j)
        MathTools::daxpy(res, -alpha, Ap);
        //z(j+1) := inv(M) * r(j+1)
        M.apply(z, res);

        double curr_resz = MathTools::innerProd(z, res);
        //beta(j) := tr(r(j+1))*z(j+1) / (tr(r(j)) * z(j))
        double beta = curr_resz/resz;

        //p(j+1) := z(j+1) + beta(j) * p(j)
        searchDirn = MathTools::vectorAdd(z, MathTools::scalarMultiply(beta, searchDirn));

        r = MathTools::L2Norm(res);

        resz = curr_resz;

        iter++;
    }

    return iter;
}

int LinearSolvers::PBiCGStab
(
    const Matrix& A,
    const Preconditioners& M,
    std::vector<double>& x,
    const std::vector<double>& rhs,
    double tolerance
)
{
    int maxIter = 10000, iter = 0;

    int n = x.size();
    vector<double> p_(n, 0.0), v(n, 0.0), s(n, 0.0), s_(n, 0.0), t(n, 0.0);

    //brief PBiCGStab (Barrett R. - Templates for solution of Linear Systems,
    //                 Building block for iterative methods)

    //First iteration is performed out of loop
    //r(0) = b - A*x(0)
    vector<double> res = MathTools::vectorSub(rhs, A.SpMV(x));
    //choosing r~ = r(0)
    vector<double> res_ = res;

    //p1 = r0 (initial search direction taken to be r0)
    vector<double> searchDirn = res;

    //rho0 = r~^T * r(0) 
    double currRho = MathTools::innerProd(res_, res);

    //inital residual
    double r0 = MathTools::L2Norm(res), r = r0;
    if(r0 < std::numeric_limits<double>::epsilon()) return iter;

    while(iter < maxIter){

        //preconditioning search direction
        M.apply(p_, searchDirn);

        //v = A*p (storing just for convinience)
        v = A.SpMV(p_);

        //step length alpha = r~^T * r / r~^T * A*p = rho / r~^T * v
        double alpha = currRho/MathTools::innerProd(res_, v);

        //s = r - alpha *  v
        s = res;
        MathTools::daxpy(s, -alpha, v);

        //checking the norm of s 
        if(MathTools::L2Norm(s)/r0 < tolerance){
            //if it's small enough, updating x and exiting
            MathTools::daxpy(x, alpha, p_);
            return iter;
        }

        //preconditioning s
        M.apply(s_, s);

        //t = A * preconditioned s
        t = A.SpMV(s_);

        //omega = t^T * s / t^T * t
        double omega = MathTools::innerProd(t, s)/MathTools::innerProd(t, t);

        if(abs(omega) < 1e-15){//for continuation it is necessary omega =/ 0
            cout << "PBiCGStab method fails\n";
            return iter;
        }

        //updating x
        //x = x + alpha * p_ + omega * s_
        MathTools::daxpy(x, alpha, p_);
        MathTools::daxpy(x, omega, s_);

        //updating the residual 
        //r = s - omega * t
        res = s;
        MathTools::daxpy(res, -omega, t);

        //calculate the Norm
        r = MathTools::L2Norm(res);

        if(r/r0 < tolerance) return iter;//no need to update next direction if converged

        //to update the search direction for next iteration
        double prevRho = currRho;
        currRho = MathTools::innerProd(res_, res);

        if(abs(currRho) < 1e-15){
            cout << "PBiCGStab Method failed\n";
            return iter;
        }

        //beta = rho(i)/rho(i-1) * (alpha/omega)
        double beta = (currRho/prevRho) * (alpha/omega);

        //p = p - omega * v
        MathTools::daxpy(searchDirn, -omega, v);
        //p = res + beta * p
        searchDirn = MathTools::vectorAdd(res, MathTools::scalarMultiply(beta, searchDirn));
        
        iter++;
    }

    return iter;
}

int LinearSolvers::steepestDescent
(
    const Matrix& A,
    vector<double>& x,
    const vector<double>& rhs, 
    double tolerance
)
{
    int maxIter = 10000, iter = 0;

    vector<double> Ar(x.size());
    vector<double> res = MathTools::vectorSub(rhs, A.SpMV(x));

    double r0 = MathTools::L2Norm(res), r = r0;
    if(r0 < std::numeric_limits<double>::epsilon()) return iter;

    while(iter < maxIter && r/r0 > tolerance){
        Ar = A.SpMV(res);
        double alpha = MathTools::innerProd(res, res);
        alpha /= MathTools::innerProd(res, Ar);

        //update x
        MathTools::daxpy(x, -alpha, res);

        //update residual
        MathTools::daxpy(res, alpha, Ar);

        //calcualte norm
        r = MathTools::L2Norm(res);

        iter++;
    }

    return iter;
}