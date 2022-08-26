#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
using namespace std;

class matrix
{
private:
public:
    int num_rows;
    int num_cols;
    double **M; 
    matrix();
    matrix(int, int);
    void print();
    int get_m();
    int get_n();

    matrix operator+(matrix);
    matrix operator-(matrix);
    matrix operator*(matrix);

    matrix operator+(double);
    matrix operator-(double);
    matrix operator*(double);
    matrix operator/(double);

    // friend matrix operator*(double, matrix);

    // friend matrix operator+(matrix, double);
    // friend matrix operator-(matrix, double);
    // friend matrix operator*(matrix, double);
    friend vector<double> operator*(matrix, vector<double>);

    matrix  transpose();
    void set_all(double);

    
};

matrix::matrix()
{
    num_rows = 0;
    num_cols = 0;
}

matrix::matrix(int m, int n)
{
    num_rows = m;
    num_cols = n;
    M = new double *[num_rows];
    for (int i = 0; i < num_rows; i++)
    {
        M[i] = new double[num_cols];
    }

    for (int i = 0; i < num_rows; i++)
    {
        for (int j = 0; j < num_cols; j++)
        {
            M[i][j] = 0;
        }
    }
}

int matrix ::get_m()
{
    return this->num_rows;
}

int matrix ::get_n()
{
    return this->num_cols;
}

void matrix::print()
{
    for (int i = 0; i < num_rows; i++)
    {
        for (int j = 0; j < num_cols; j++)
        {   
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// ----------------------  <<<< Operator Overloading  Matrix * vector >>> -------------------//
vector<double> operator*(matrix M, vector<double> u)
{
    vector<double> v(u.size());

    if(M.get_n() == u.size())
    {
        for(int i = 0; i < M.get_m(); i++)
        {
            for(int j = 0; j < M.get_n(); j++)
            {
                v[i] += M.M[i][j]*u[j];
            }
        }
    }
    else
    {
        cout << "invalid matrix vector multiplication, check dimensions." << endl;
    }

}

// ----------------------  <<<<  (Matrix,Matrix) Operator Overloading >>> -------------------//

matrix matrix::operator+(matrix A)
{
    matrix P;
    if(num_rows != A.num_rows || num_cols != A.num_cols)
    {
        cout << "Matrix addition is invalid\n";
        exit(0);
    }
    else
    {
        P = matrix(num_rows, num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<num_cols;j++)
                P.M[i][j] = M[i][j] + A.M[i][j];
        return P;
    }
}

matrix matrix::operator-(matrix A)
{
    matrix P;
    if(num_rows != A.num_rows || num_cols != A.num_cols)
    {
        cout << "Matrix subtraction is invalid\n";
        exit(0);
    }
    else
    {
        P = matrix(num_rows, num_cols);
        for(int i=0;i<num_rows;i++)
        {
            for(int j=0;j<num_cols;j++)
            {   
                P.M[i][j] = M[i][j] - A.M[i][j];
            }
        }    
        return P;
    }
}


matrix matrix::operator*(matrix A)
{
    matrix P;
    if (this->num_cols != A.num_rows)
    {
        cout << "Invalid matrix multiplication ! \n";
        exit(0);
    }
    else
    {
        P = matrix(this->num_rows, A.num_cols);
        for(int i=0;i<this->num_rows;i++)
            for(int j=0;j<A.num_cols;j++)
            {
                P.M[i][j]=0.0;
                for(int k=0;k<this->num_cols;k++)
                {
                    P.M[i][j] = P.M[i][j] + this->M[i][k]*A.M[k][j];
                }
            }
            cout <<endl;
        return P;
    }
}

// ----------------------  <<<<  (Matrix,double) Operator Overloading >>> -------------------//

matrix matrix::operator+(double a)
{
    matrix P;
    P = matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] + a;
    return P;
    
}

matrix matrix::operator-(double a)
{
    matrix P;
    P = matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] - a;
    return P;
    
}
matrix matrix::operator*(double a)
{
    matrix P;
    P = matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] * a;
    return P;
    
}
matrix matrix::operator/(double a)
{
    if(a!=0)
    {
    matrix P;
    P = matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] / a;
    return P;
    }
    else
    {
        cout << "Cannot perform division by zero"<< endl;
    }
    
}


// ------------------------------

matrix matrix :: transpose()
{
    matrix A(num_cols, num_rows);

    for(int i = 0; i < A.get_m(); i++)
    {
        for(int j = 0; j < A.get_n(); j++)
        {
            A.M[i][j] = M[j][i];
        }
    }
    return A;
}

void matrix ::set_all(double value)
{
    for(int i = 0; i < num_rows; i ++)
    {
        for(int j = 0; j < num_cols; j++)
        {
            M[i][j] = value;
        }
    }
}

#endif