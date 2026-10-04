#include "MathTools.h"

#include <cmath>

using namespace std;

double MathTools::L2Norm
(
    const vector<double>& res
)
{
    int n = res.size();
    double l2 = 0;
    for(int i = 0; i < n; ++i){
        l2 += res[i] * res[i];
    }

    l2 = sqrt(l2);

    return l2;
}

vector<double> MathTools::vectorSub
(
    const vector<double>& v1,
    const vector<double>& v2
)
{
    vector<double> r;
    int n = v1.size();
    r.assign(n, 0.0);

    for(int i = 0; i<n; ++i){
        r[i] = v1[i] - v2[i];
    }

    return r;
}

vector<double> MathTools::vectorAdd
(
    const vector<double>& v1,
    const vector<double>& v2
)
{
    vector<double> r;
    int n = v1.size();
    r.assign(n, 0.0);

    for(int i = 0; i<n; ++i){
        r[i] = v1[i] + v2[i];
    }

    return r;
}

double MathTools::innerProd
(
    const vector<double>& v1,
    const vector<double>& v2
)
{
    double prod = 0.0;

    for(int i = 0; i < v1.size(); ++i){
        prod += v1[i] * v2[i];
    }

    return prod;
}

vector<double> MathTools::scalarMultiply
(
    double scalar,
    const vector<double>& v
)
{
    vector<double> sv(v.size());

    for(int i = 0; i < v.size(); ++i){
        sv[i] = scalar * v[i];
    }

    return sv;
}

void MathTools::daxpy
(
    std::vector<double>& x,
    const double a,
    const std::vector<double>& y
)
{
    for(int i = 0; i < x.size(); ++i){
        x[i] = x[i] + a*y[i];
    }
}