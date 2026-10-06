#include "Smoothers.h"

using namespace std;

void Smoothers::gaussSeidel
(
    Matrix& A,
    const vector<double>& rhs,
    vector<double>& x,
    int maxIter
)
{
    auto rowPtr = A.getrowPtr();
    auto col = A.getcol();
    auto values = A.getvalues();
    
    int n = rowPtr.size() - 1;
    vector<double> invDiag(n, 0.0);

    for(int i = 0; i < n; ++i){
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            if(i == col[j]){
                invDiag[i] = 1/values[j];
                break;
            }
        }
    }

    for(int iter = 0; iter < maxIter; ++iter){
        for(int i = 0; i < n; ++i){
            double sum = rhs[i];
            for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
                if(i == col[j]) continue;

                sum -= values[j] * x[col[j]];
            }

            x[i] = sum * invDiag[i];
        }
    }
}