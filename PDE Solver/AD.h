#ifndef AD_H
#define AD_H

#include <iostream>
#include <string>
#include <vector>
#include "Matrix.h"
using namespace std;
#include <cmath>

//------------------------------------------------------    AD class    -------------------------------------------------//

class AD
{
private:
    double f;
    vector<double> df;
    int id;

public:
    AD();
    AD(double, int);
    AD(double, int, int);
    void setIndVar(int size);
    double getf();
    vector<double> getdf();
    double getDf(int);
    int getID();
    vector<double> getGradient();
    friend matrix getJacobian(vector<AD>, int n);

    // ---------------------------  Operator Overloading on f and g functions, c is a constant ---------------------------//

    // -------------------------  Operator Overloading for  f+g, f-g, f*g, f/g, f^c, declaration -------------------------//

    AD operator+(AD);  // f + g
    AD operator-(AD);  // f - g
    AD operator*(AD);  // f * g
    AD operator/(AD);  // f / g
    AD operator^(int); // f^c

    // --------------------------  Operator Overloading for  c+f, c-f, c*f, c/f, c^f, declaration ------------------------//

    friend AD operator+(double, AD); // c + f
    friend AD operator-(double, AD); // c - f
    friend AD operator*(double, AD); // c * f
    friend AD operator/(double, AD); // c / f
    friend AD operator^(double, AD); // c ^ f

    // ---------------------------  Operator Overloading for  f+c, f-c, f*c, f/c, f^g, declaration -----------------------//

    AD operator+(double); // f + c
    AD operator-(double); // f - c
    AD operator*(double); // f * c
    AD operator/(double); // f / c
    AD operator^(AD);     // f ^ g

    /*-------------------- Declaring friend functions for the different mathematical fuctions ----------------------------

									sin(f), cos(f), tan(f), cosec(f), sec(f), cot(f).
									arcsin(f), arccos(f), arctan(f), sinh(f), cosh(f), tanh(f).
									log(f), exp(f), abs(f) 											 						*/

    friend AD sin(AD);    // sin()
    friend AD cos(AD);    // cos()
    friend AD tan(AD);    // tan()
    friend AD cosec(AD);  // cosec()
    friend AD sec(AD);    // sec()
    friend AD cot(AD);    // cot()
    friend AD arcsin(AD); // arcsin()
    friend AD arccos(AD); // arccos()
    friend AD arctan(AD); // arctan()
    friend AD sinh(AD);   // sinh()
    friend AD cosh(AD);   // cosh()
    friend AD tanh(AD);   // tanh()
    friend AD log(AD);    // log()
    friend AD exp(AD);    // exp()
    friend AD abs(AD);    // abs()
};

// ----------------------------------- Class AD Constructors and Member fuctions ----------------------------------------//

AD::AD()
{   // Default AD constructor
    f = 0;
}

AD::AD(double value, int size)
{   // AD constructor for converting a double constant value into an object of AD type
    this->f = value;
    this->df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        this->df[i] = 0;
    }
}

AD::AD(double value, int ID, int size)
{   // AD constructor to make an independent variable,
    // just give the value, ID, size of the ind var, size is the size of df vector<double>
    this->f = value;
    this->id = ID;
    this->df = vector<double>(size);
    for (int i = 0; i < this->id; i++)
    {
        this->df[i] = 0;
    }
    this->df[this->id] = 1;
    for (int i = this->id + 1; i < size; i++)
    {
        this->df[i] = 0;
    }
}

void AD ::setIndVar(int size)
{   // sets the independent variables, creates a vector<double> of derivatives, and set the derivative wrt the same variable as 1
    // otherwise set 0
    this->df = vector<double>(size);
    for (int i = 0; i < this->id; i++)
    {
        this->df[i] = 0;
    }
    this->df[this->id] = 1;
    for (int i = this->id + 1; i < size; i++)
    {
        this->df[i] = 0;
    }
}

double AD::getf()
{ // returns the value of function f
    return this->f;
}

vector<double> AD::getdf()
{ // returns the value of function f
    return this->df;
}

double AD::getDf(int index)
{ // gives the partial derivative of function wrt a variable(based on index)  for index = 0, df/dx; index = 1, df/dy;
    return df[index];
}

int AD::getID()
{ // gives the partial derivative of function wrt a variable(based on index)  for index = 0, df/dx; index = 1, df/dy;
    return this->id;
}

//------------------------------------------    getGradient function     ---------------------------------------//
vector<double> AD::getGradient()
{   
    int size = this->df.size();
    vector<double> gradient(size);
    for (int i = 0; i < size; i++)
    {
        gradient[i] = this->df[i];
        // cout << " gradient wrt arr[i] element :"<< this->df[i] << endl;
    }
    return gradient;
}

// -------------------------------     getJacobian function  friend of AD class     ----------------------------//
// Pls check this, could create problem later on
matrix getJacobian(vector<AD> funList)
{   // n is the number of elements in the funList array.
    int n = funList.size();
    matrix M = matrix(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            M.M[i][j] = funList[i].getDf(j); // look at the index of Df, it should match with the id of AD object
            // cout << "derivative of f" << i << " wrt " << " arr["<<j <<"] element : " <<M.getElement(i,j) <<endl;
        }
    }
    return M;
}

//-----------------------------------  Operator Overloading for  f+g, f-g, f*g, f/g, f^c, Definitions  ---------------------------------//

AD AD::operator+(AD g)
{   // f + g
    AD h;
    h.f = this->f + g.f;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = this->df[i] + g.df[i];
    }
    return h;
}

AD AD::operator-(AD g)
{   // f - g
    AD h;
    h.f = this->f - g.f;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = this->df[i] - g.df[i];
    }
    return h;
}

AD AD::operator*(AD g)
{   // f * g
    AD h;
    h.f = this->f * g.f;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = this->f * g.df[i] + g.f * this->df[i];
    }
    return h;
}

AD AD::operator/(AD g)
{   // f / g
    AD h;
    h.f = this->f / g.f;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = (this->df[i] * g.f - g.df[i] * this->f) / (g.f * g.f);
    }
    return h;
}

AD AD::operator^(int c)
{   // f ^ c
    AD h;
    h.f = pow(this->f, c);
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = c * pow(this->f, c - 1);
    }
    return h;
}

// -----------------------------------  Operator Overloading for  c+f, c-f, c*f, c/f, c^f, Definitions  --------------------------------//

AD operator+(double c, AD g)
{   // c + f
    AD h;
    h.f = c + g.f;
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = 0 + g.df[i];
    }
    return h;
}

AD operator-(double c, AD g)
{   // c - f
    AD h;
    h.f = c - g.f;
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = 0 - g.df[i];
    }
    return h;
}

AD operator*(double c, AD g)
{   // c * f
    AD h;
    h.f = c * g.f;
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = c * g.df[i] + 0 * g.f;
    }
    return h;
}

AD operator/(double c, AD g)
{   // c / f
    AD h;
    h.f = c / g.f;
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = (0 * g.f - g.df[i] * c) / (g.f * g.f);
    }
    return h;
}

AD operator^(double c, AD g)
{   // c ^ f
    AD h;
    h.f = pow(c, g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = pow(c, g.f) * log(c) * g.df[i];
    }
    return h;
}

// ----------------------------------  Operator Overloading for  f+c, f-c, f*c, f/c, f^g, Definitions  ---------------------------------//

AD AD::operator+(double c)
{   // f + c
    AD h;
    h.f = this->f + c;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = this->df[i] + c;
    }
    return h;
}

AD AD::operator-(double c)
{   // f - c
    AD h;
    h.f = this->f - c;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = this->df[i] - c;
    }
    return h;
}

AD AD::operator*(double c)
{   // f * c
    AD h;
    h.f = this->f * c;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = this->f * 0 + this->df[i] * c;
    }
    return h;
}

AD AD::operator/(double c)
{   // f / c
    AD h;
    h.f = this->f / c;
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = (this->df[i] * c - 0 * this->f) / (c * c);
    }
    return h;
}

AD AD::operator^(AD g)
{   // f ^ g
    AD h;
    h.f = pow(this->f, g.f);
    int size = this->df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = pow(this->f, g.f) * (g.f * (this->df[i] / this->f) + g.df[i] * log(abs(this->f)));
    }
    return h;
}

//--------------------------- Defining friend functions for the different mathematical fuctions ---------------------------//

AD sin(AD g)
{   // sin() 
    AD h;
    h.f = sin(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = cos(g.f) * g.df[i];
    }
    return h;
}

AD cos(AD g)
{   // cos() 
    AD h;
    h.f = cos(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = -1 * sin(g.f) * g.df[i];
    }
    return h;
}

AD tan(AD g)
{   // tan() 
    AD h;
    h.f = tan(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = (1 / cos(g.f)) * (1 / cos(g.f)) * g.df[i];
    }
    return h;
}

AD cosec(AD g)
{   // cosec() 
    AD h;
    h.f = (1 / sin(g.f));
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = -1 * (1 / sin(g.f)) * (1 / tan(g.f)) * g.df[i];
    }
    return h;
}

AD sec(AD g)
{   // sec() 
    AD h;
    h.f = (1 / cos(g.f));
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = (1 / cos(g.f)) * tan(g.f) * g.df[i];
    }
    return h;
}

AD cot(AD g)
{   // cot() 
    AD h;
    h.f = (1 / tan(g.f));
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = -1 * (1 / sin(g.f)) * (1 / sin(g.f)) * g.df[i];
    }
    return h;
}

AD arcsin(AD g)
{   // arcsin() 
    AD h;
    h.f = asin(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        h.df[i] = (1 / sqrt(1 - g.f * g.f)) * g.df[i];
    }
    return h;
}

AD arccos(AD g)
{   // arccos() 
    AD h;
    h.f = acos(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = -1 * (1 / sqrt(1 - g.f * g.f)) * g.df[i];
    }
    return h;
}

AD arctan(AD g)
{   // arctan() 
    AD h;
    h.f = atan(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = (1 / (1 + g.f * g.f)) * g.df[i];
    }
    return h;
}

AD sinh(AD g)
{   // sinh() 
    AD h;
    h.f = sinh(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = cosh(g.f) * g.df[i];
    }
    return h;
}

AD cosh(AD g)
{   // cosh() 
    AD h;
    h.f = cosh(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = sinh(g.f) * g.df[i];
    }
    return h;
}

AD tanh(AD g)
{   // tanh() 
    AD h;
    h.f = tanh(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = (1 - (tanh(g.f) * tanh(g.f))) * g.df[i];
    }
    return h;
}

AD log(AD g)
{   // log() 
    AD h;
    h.f = log(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = (1 / g.f) * g.df[i];
    }
    return h;
}

AD exp(AD g)
{   // exp() 
    AD h;
    h.f = exp(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = exp(g.f) * g.df[i];
    }
    return h;
}

AD abs(AD g)
{   // abs() 
    AD h;
    h.f = abs(g.f);
    int size = g.df.size();
    h.df = vector<double>(size);
    for (int i = 0; i < size; i++)
    {
    	h.df[i] = (g.f / abs(g.f)) * g.df[i];
    }
    return h;
}

#endif