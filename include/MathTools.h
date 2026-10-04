#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <vector>

class MathTools{
public:

    static double L2Norm
    (
        const std::vector<double>& res
    );

    static std::vector<double> vectorSub
    (
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );

    static std::vector<double> vectorAdd
    (
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );

    static double innerProd
    (
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );

    static std::vector<double> scalarMultiply
    (
        double scalar,
        const std::vector<double>& v
    );

    static void daxpy
    (
        std::vector<double>& x,
        const double a,
        const std::vector<double>& y
    );

};
#endif