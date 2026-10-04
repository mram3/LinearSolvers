#include "Matrix.h"

#include <algorithm>

using namespace std;

void Matrix::addCoeff
(
    int row,
    int colm,
    double value
)
{
    rowsCOO.push_back(row);
    col.push_back(colm);
    values.push_back(value);
}

struct nodes{
    int col;
    double val;

    bool operator <(const nodes& other) const{
        return col < other.col;
    }
};

void Matrix::compressToCSR()
{
    int colCount = 0, count = 0;

    rowPtr.push_back(count);
    for(int i = 0; i < rowsCOO.size(); ++i){

        if(i == rowsCOO.size()-1 || rowsCOO[i]!=rowsCOO[i+1]){
            int num_elements = colCount + 1;

            vector<nodes> row_data(num_elements);
            for(int j = 0; j < num_elements; ++j){
                row_data[j] = {col[count+j], values[count+j]};
            }

            sort(row_data.begin(), row_data.end());
            
            for(int j = 0; j < num_elements; ++j){
                col[count+j] = row_data[j].col;
                values[count+j] = row_data[j].val;
            }

            count += num_elements;
            rowPtr.push_back(count);
            colCount = 0;
        }
        else{
            colCount++;
        }
    }

    rowsCOO.clear();
    rowsCOO.shrink_to_fit();
}

vector<double> Matrix::SpMV
(
    const vector<double>& x
) const
{
    int n = rowPtr.size() - 1;
    vector<double> r(n);

    for(int i = 0; i<n; ++i){
        double sum = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            sum += values[j] * x[col[j]];
        }
        r[i] = sum;
    }

    return r;
}

const vector<int>& Matrix::getrowPtr() const
{
    return rowPtr;
}

const vector<int>& Matrix::getcol() const
{
    return col;
}

const vector<double>& Matrix::getvalues() const
{
    return values;
}