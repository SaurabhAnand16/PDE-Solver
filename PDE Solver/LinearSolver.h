#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <iostream>
using namespace std;
#include <cmath>
#include <vector>
#include "Matrix.h"

// 1) Gauss Elimination Method
vector<double> GaussElimination_Solver(matrix M, vector<double> u)
{   
    vector<double> x(u.size());
    double temp;
    int i, j, k;
    for (i = 0; i < u.size(); i++)
    {
        temp = fabs(M.M[i][i]);
        k = i;
        for (j = i + 1; j < u.size(); j++)
            if (temp < fabs(M.M[j][i]))
            {
                temp = fabs(M.M[j][i]);
                k = j;
            }
        if (fabs(M.M[k][i]) < 0.00001)
        {
            cout << "The matrix is singular: The system has either no solution or infinitely many solution";
            exit(0);
        }
        if (k != i)
        {
            for (j = 0; j < u.size(); j++)
            {
                temp = M.M[k][j];
                M.M[k][j] = M.M[i][j];
                M.M[i][j] = temp;
            }
            temp = u[k];
            u[k] = u[i];
            u[i] = temp;
        }
        for (j = i + 1; j < u.size(); j++)
        {
            temp = M.M[j][i] / M.M[i][i];
            for (k = 0; k < u.size(); k++)
                M.M[j][k] = M.M[j][k] - temp * M.M[i][k];
            u[j] = u[j] - temp * u[i];
        }
    }
    x[u.size() - 1] = u[u.size() - 1] / M.M[u.size() - 1][u.size() - 1];
    for (i = u.size() - 2; i >= 0; i--)
    {
        x[i] = u[i];
        for (j = i + 1; j < u.size(); j++)
            x[i] = x[i] - M.M[i][j] * x[j];
        x[i] = x[i] / M.M[i][i];
    }
    return x;
}

// 2) LU Decomposition Method
vector<double> LU_Decomposition_Solver(matrix M, vector<double> u)
{
    vector<double> x(u.size()), z(u.size());
    matrix l(u.size(), u.size()), u1(u.size(),u.size());
    int i,j,k;
    
    for(i=0;i<u.size();i++)
        l.M[i][0] = M.M[i][0];
    for(j=1;j<u.size();j++)
        u1.M[0][j]= M.M[0][j]/l.M[0][0];
    for(i=0;i<u.size();i++)
        u1.M[i][i] = 1 ;
    for(i=1;i<u.size();i++)
        for(j=1;j<u.size();j++)
            if(i>=j)
            {
                l.M[i][j] =M.M[i][j];
                for(k=0;k<=j-1;k++){
                    
                    l.M[i][j] = l.M[i][j]-l.M[i][k]*u1.M[k][j];
                }
            }
            else
            {
                u1.M[i][j] =M.M[i][j];
                for(k=0;k<=i-1;k++)
                    u1.M[i][j]= u1.M[i][j]-l.M[i][k]*u1.M[k][j];
                u1.M[i][j]= u1.M[i][j]/l.M[i][i];
            }
    cout<<"The lower triangular matrix L:"<<endl;
    for(i=0;i<u.size();i++)
    {
        for(j=0;j<=i;j++)
            cout<<"\t"<<l.M[i][j];
        cout<<endl;
    }
    cout<<"\nThe upper triangular matrix U:"<<endl;
    for(i=0;i<u.size();i++)
    {
        for(j=0;j<i;j++)
            cout<<"\t";
        for(j=i;j<u.size();j++)
            cout<<"\t"<<u1.M[i][j];
        cout<<endl;
    }
    z[0]=u[0]/l.M[0][0];
    for(i=1;i<u.size();i++)
    {
        z[i]=u[i];
        for(j=0;j<=i-1;j++)
            z[i]-=l.M[i][j]*z[j];
        z[i]/=l.M[i][i];
    }
    x[u.size()-1]=z[u.size()-1];
    for(i=u.size()-2;i>=0;i--)
    {
        x[i]=z[i];
        for(j=i+1;j<u.size();j++)
            x[i]-=u1.M[i][j]*x[j];
    }
    return x;
}

// 3) TriDiagonal Method
vector<double> TriDiagonal_Solver(matrix M, vector<double> u)
{
    vector<double> x(u.size());
    double A[u.size()-1],B[u.size()],C[u.size()-1],D[u.size()];
    double C_star[u.size()-1],D_star[u.size()];
    int i,j;
    
    for(i=0;i<u.size()-2;i++)
        for(j=i;j<u.size()-2;j++)
    {
        if(M.M[i][j+2]!=0)
        {
            cout<<"Method can't be applied";
            exit(0);
        }

    }
    for(i=u.size()-1;i>1;i--)
        for(j=i;j>1;j--)
    {
        if(M.M[i][j-2]!=0)
        {
            cout<<"Method can't be applied";
            exit(0);
        }
    }
    for(i=0;i<u.size();i++)
    {
        B[i]=M.M[i][i];
        D[i]=u[i];
    }
    for(i=1;i<u.size();i++)
    {
       A[i]=M.M[i][i-1];
    }
    for(i=0;i<u.size()-1;i++)
    {
       C[i]=M.M[i][i+1];
    }

    C_star[0]=C[0]/B[0];
    D_star[0]=D[0]/B[0];
    for(i=1;i<u.size();i++)
    {
        C_star[i]=C[i]/(B[i]-A[i]*C_star[i-1]);
        D_star[i]=(D[i]-A[i]*D_star[i-1])/(B[i]-A[i]*C_star[i-1]);
    }
    x[u.size()-1]=D_star[u.size()-1];
    for(i=u.size()-2;i>=0;i--)
    {
        x[i]=D_star[i]-C_star[i]*x[i+1];
    }
      return x;
}

// 4) Gauss Jacobi Method

vector<double> Gauss_Jacobi_Solver(matrix M, vector<double> u)
{
    vector<double> x,xn;
    int i,j,flag;
    double sum,eps;
    
    x=vector<double>(u.size());
    xn=vector<double>(u.size());
    cout<<"Enter the accuracy u want\n";
    cin>>eps;
    for(i=0;i<u.size();i++)
        x[i]=0;
    flag=0;
    for(i=0;i<u.size();i++)
    {
        sum=0;
        for(j=0;j<u.size();j++)
            if(i!=j)
                sum+=fabs(M.M[i][j]);
        if(sum>fabs(M.M[i][i]))
            flag=1;
    }
    if(flag==1)
    {
        flag=0;
        for(j=0;j<u.size();j++)
        {
            sum=0;
            for(i=0;i<u.size();i++)
                if(i!=j)
                    sum+=fabs(M.M[i][j]);
            if(sum>fabs(M.M[j][j]))
                flag=1;
        }
    }
    if(flag==1)
    {
        cout<<"The co-efficient matrix is not diagonally dominant\n";
        cout<<"The Gauss-Jacobi method doesn't converge surely\n";
        exit(0);
    }
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i];
            for(j=0;j<u.size();j++)
                if(j!=i)
                    sum-=M.M[i][j]*x[j];
            xn[i]=sum/M.M[i][i];
        }
        flag=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                flag=1;
        if(flag==1)
            for(i=0;i<u.size();i++)
                x[i]=xn[i];
    }while(flag==1);
    return xn;
}
// 5) Gauss Seidal Method
vector<double> Gauss_Seidal_Solver(matrix M, vector<double> u)
{
    vector<double> x,xn;
    int test,i,j;
    double g,sum,eps=0.0001;
    //cout<<"SATYA";
    for(i=0;i<u.size();i++){
       for(j=0;j<u.size();j++){
         cout<<M.M[i][j]<<" ";
       }
       cout<<"\n";
    }

    x=vector<double>(u.size());
    xn=vector<double>(u.size());
    for(i=0;i<u.size();i++)
    {
        x[i] = 0;
    }
    
    for(i=0;i<u.size();i++)
        x[i]=0;

        
    test=0;
    for(i=0;i<u.size();i++)
    {
        sum=0;
        for(j=0;j<u.size();j++)
            if(i!=j)
                sum+=fabs(M.M[i][j]);
        if(sum>fabs(M.M[i][i]))
            test=1;
    }
 // Diagonally dominance verification by columns.
    if(test==1)
    {
        test=0;
        for(j=0;j<u.size();j++)
        {
            sum=0;
            for(i=0;i<u.size();i++)
                if(i!=j)
                    sum+=fabs(M.M[i][j]);
            if(sum>fabs(M.M[j][j]))
                test=1;
        }
    }

    if(test==1)
    {
        cout<<"The co-efficient matrix is not diagonally dominant\n";
        cout<<"The Gauss-Seidel method doesn't converge surely\n";
        exit(0);
    }
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i];
            for(j=0;j<u.size();j++)
            {
                if(j<i)
                    sum-=M.M[i][j]*xn[j];
                else if(j>i)
                    sum-=M.M[i][j]*x[j];
            }
            xn[i]=sum/M.M[i][i];
        }
        test=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                test=1;
        if(test==1)
            for(i=0;i<u.size();i++)
                x[i] = xn[i];
    }while(test==1);
    return xn;
}

// 6) Successive Over Relaxation Method

vector<double> SOR_Solver(matrix M, vector<double> u){
    vector<double> x,xn;
    int i,j,flag;
    double sum,eps=0.001,w=1;
    
    x=vector<double>(u.size());
    xn=vector<double>(u.size());
    
    
    for(i=0;i<u.size();i++)
        x[i]=0;
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i]*w+M.M[i][i]*x[i];
            for(j=0;j<u.size();j++)
            {
                if(j<i)
                    sum-=M.M[i][j]*xn[j]*w;
                else if(j>=i)
                    sum-=M.M[i][j]*x[j]*w;
                xn[i]=sum/M.M[i][i];
            }
        }
        flag=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                flag=1;
        if(flag==1)
            for(i=0;i<u.size();i++)
                x[i]=xn[i];
    }while(flag==1);
    return xn;
}





#endif