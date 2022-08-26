#ifndef NONLINEAR_SOLVER_H
#define NONLINEAR_SOLVER_H

#include <iostream>
using namespace std;
// #include "infix_to_postfix.h"
#include "Discretizer.h"
#include "Matrix.h"
#include "LinearSolver.h"

class PDE
{
private:
    string A, B, C, D, E;
    stringList postfix_A, postfix_B, postfix_C, postfix_D, postfix_E;

public:
    Discretize domain;
    PDE();
    PDE(string, string, string, string, string, Discretize);
    vector<AD> F;
    vector<AD> evaluate_PDE(vector<double>);
};

PDE ::PDE()
{
    this->A = "";
    this->B = "";
    this->C = "";
    this->D = "";
    this->E = "";
}

PDE ::PDE(string expression_A, string expression_B, string expression_C, string expression_D, string expression_E, Discretize domain)
{
    this->A = expression_A;
    this->B = expression_B;
    this->C = expression_C;
    this->D = expression_D;
    this->E = expression_E;
    this->domain = domain;
    this->F = vector<AD>(domain.ind_var_Points.size());

    stringList infix_A, infix_B, infix_C, infix_D, infix_E;
    this->postfix_A = infix_A.infix_to_postfix(A);
    this->postfix_B = infix_B.infix_to_postfix(B);
    this->postfix_C = infix_C.infix_to_postfix(C);
    this->postfix_D = infix_D.infix_to_postfix(D);
    this->postfix_E = infix_E.infix_to_postfix(E);
}

vector<AD> PDE ::evaluate_PDE(vector<double> U_value)
{
    // update the u_value of ind_var_point vector with the new U_value(taken as input to the this function)
    for (int k = 0; k < domain.ind_var_Points.size(); k++)
    {
        AD new_u(U_value[k], k, U_value.size());
        domain.ind_var_Points[k].set_u(new_u);
    }

    for (int k = 0; k < domain.ind_var_Points.size(); k++)
    {
        int i = domain.ind_var_Points[k].get_i();
        int j = domain.ind_var_Points[k].get_j();

        double x_val = domain.ind_var_Points[k].get_x_val();
        double y_val = domain.ind_var_Points[k].get_y_val();

        int no_of_ind_vars = domain.ind_var_Points.size();

        AD u_val = domain.ind_var_Points[k].get_u();

        AD coeff_A = postfix_A.evaluate(x_val, y_val, no_of_ind_vars, u_val);
        AD coeff_B = postfix_B.evaluate(x_val, y_val, no_of_ind_vars, u_val);
        AD coeff_C = postfix_C.evaluate(x_val, y_val, no_of_ind_vars, u_val);
        AD coeff_D = postfix_D.evaluate(x_val, y_val, no_of_ind_vars, u_val);
        AD coeff_E = postfix_E.evaluate(x_val, y_val, no_of_ind_vars, u_val);

        AD d2u_by_dx2 = (domain.get_Point(i + 1, j).get_u() - 2 * domain.get_Point(i, j).get_u() + domain.get_Point(i - 1, j).get_u()) / (domain.delta_x * domain.delta_x);
        AD d2u_by_dy2 = (domain.get_Point(i, j + 1).get_u() - 2 * domain.get_Point(i, j).get_u() + domain.get_Point(i, j - 1).get_u()) / (domain.delta_y * domain.delta_y);
        AD du_by_dx = (domain.get_Point(i + 1, j).get_u() - domain.get_Point(i - 1, j).get_u()) / (2 * domain.delta_x);
        AD du_by_dy = (domain.get_Point(i, j + 1).get_u() - domain.get_Point(i, j - 1).get_u()) / (2 * domain.delta_y);

        AD f_k = coeff_A * d2u_by_dx2 + coeff_B * d2u_by_dy2 + coeff_C * du_by_dx + coeff_D * du_by_dy + coeff_E;

        this->F[k] = f_k;
    }
    return F;
}

Discretize Newton_Solver(PDE gen_PDE, int n);
Discretize Broyden_Solver(PDE gen_PDE, int n);

Discretize Solve(PDE gen_PDE, int NonLinearChoice, int LinearChoice)
{
    Discretize solution;
    switch (NonLinearChoice)
    {
    case 1:
        solution = Newton_Solver(gen_PDE, LinearChoice);
        break;
    case 2:
        solution = Broyden_Solver(gen_PDE, LinearChoice);
        break;
    default:
        cout << "Please choose an appropriate Non-Linear Solver Method." << endl;
        break;
    }
    return solution;
}
//----------------------------------------------------------------------------------------


vector<double> makeDoubleVector(vector<AD> F_X_AD)
{
    vector<double> F_X_dobule;
    for (int i = 0; i < F_X_AD.size(); i++)
    {
        F_X_dobule.push_back(F_X_AD[i].getf());
    }
    return F_X_dobule;
}

double norm(vector<double> X_k0, vector<double> X_k1)
{
    double error_sqr = 0;
    for (int i = 0; i < X_k0.size(); i++)
    {
        error_sqr += pow((X_k1[i] - X_k0[i]), 2);
    }
    return sqrt(error_sqr);
}

matrix inverse_matrix(matrix A, int method_num)
{
    matrix inv_A(A.get_m(), A.get_n());
    int size = A.get_m();
    matrix copy_A(size, size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            copy_A.M[i][j] = A.M[i][j];
        }
    }
    if (A.get_m() != A.get_n())
    {
        cout << "Matrix is not a square matrix" << endl;
    }
    else
    {
        vector<double> x(A.get_m());
        for (int k = 0; k < A.get_m(); k++)
        {
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                     A.M[i][j] = copy_A.M[i][j] ;
                }
            }
            vector<double> b(A.get_m(), 0);
            b[k] = 1;

            switch (method_num)
            {
            case 1:
                x = GaussElimination_Solver(A, b);
                break;
            case 2:
                x = LU_Decomposition_Solver(A, b);
                break;
            case 3:
                x = TriDiagonal_Solver(A, b);
                break;
            case 4:
                x = Gauss_Jacobi_Solver(A, b);
                break;
            case 5:
                x = Gauss_Seidal_Solver(A, b);
                break;
            case 6:
                x = SOR_Solver(A, b);
                break;
            default:
                cout << "Please enter a valid number corrosponding to the Linear Solver type" << endl;
                break;
            }
            for (int j = 0; j < size; j++)
            {
                inv_A.M[j][k] = x[j];
            }
        }
    }
    return inv_A;
}

vector<double> vector_from_matrix(matrix M)
{
    int m = M.get_m();
    int n = M.get_n();
    vector<double> v((m * n), 0);
    // k = v.size();
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            v[n * i + j] = M.M[i][j];
        }
    }
    return v;
}

matrix matrix_from_vector(vector<double> v, int m, int n)
{
    matrix M(m, n);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int index = n * i + j;
            M.M[i][j] = v[n * i + j];
        }
    }
    return M;
}

//______________________________   Newton Solver fuction Definitions   ______________________________________//

Discretize Newton_Solver(PDE gen_PDE, int n)
{

    vector<double> initial_guess(gen_PDE.domain.ind_var_Points.size(), 0.1);
    vector<double> X_k0_iter = initial_guess;

    vector<AD> F_X_k0_iter;
    vector<double> F_X_k0_double;

    matrix J_X_k0_iter;
    vector<double> X_k1_iter(gen_PDE.domain.ind_var_Points.size(), 0);
    vector<double> Y_k0_iter;
    int count = 1;
    double error = 1;
    double epsilon = 0.00001;

    while (error > epsilon)
    {
        F_X_k0_iter = gen_PDE.evaluate_PDE(X_k0_iter);
        F_X_k0_double = makeDoubleVector(F_X_k0_iter);

        J_X_k0_iter = getJacobian(F_X_k0_iter);
        
        switch (n)
        {
        case 1:
            Y_k0_iter = GaussElimination_Solver(J_X_k0_iter, F_X_k0_double);
            break;
        case 2:
            Y_k0_iter = LU_Decomposition_Solver(J_X_k0_iter, F_X_k0_double);
            break;
        case 3:
            Y_k0_iter = TriDiagonal_Solver(J_X_k0_iter, F_X_k0_double);
            break;
        case 4:
            Y_k0_iter = Gauss_Jacobi_Solver(J_X_k0_iter, F_X_k0_double);
            break;
        case 5:
            Y_k0_iter = Gauss_Seidal_Solver(J_X_k0_iter, F_X_k0_double);
            break;
        case 6:
            Y_k0_iter = SOR_Solver(J_X_k0_iter, F_X_k0_double);
            break;
        default:
            cout << "Please enter a valid number corrosponding to the Linear Solver type" << endl;
            break;
        }

        for (int i = 0; i < X_k0_iter.size(); i++)
        {
            X_k1_iter[i] = X_k0_iter[i] - Y_k0_iter[i];
        }

        error = norm(X_k0_iter, X_k1_iter);

        if (error < epsilon)
        {
            break;
        }
        else
        {
            X_k0_iter = X_k1_iter;
            count++;
        }
    }
    vector<double> final_soln = X_k1_iter;

    for (int k = 0; k < gen_PDE.domain.ind_var_Points.size(); k++)
    {
        AD temp(final_soln[k], k, gen_PDE.domain.ind_var_Points.size());
        gen_PDE.domain.ind_var_Points[k].set_u(temp);
    }
    return gen_PDE.domain;
}

//____________________________________   Broyden Solver fuction Definitions   ________________________________//

Discretize Broyden_Solver(PDE gen_PDE, int n)
{
    cout << "entered Broyden solver pde" << endl;
    int no_of_elements = gen_PDE.domain.ind_var_Points.size();
    int count = 1;
    double epsilon = 0.0001;
    double error = 1;

    // step 1       
    double initial_guess = 0.1;
    vector<double> initial_guess_vec(gen_PDE.domain.ind_var_Points.size(), initial_guess);
    vector<double> X0_vec = initial_guess_vec;
    matrix X0 = matrix_from_vector(X0_vec, no_of_elements, 1);

    // step 2
    vector<AD> FX0_AD_vec = gen_PDE.evaluate_PDE(X0_vec);
    vector<double> FX0_vec = makeDoubleVector(FX0_AD_vec);
    matrix FX0 = matrix_from_vector(FX0_vec, no_of_elements, 1);

    //step 3
    matrix JX0 = getJacobian(FX0_AD_vec);
    matrix A0_inv = inverse_matrix(JX0, n);
    // step 4
    matrix X1 = X0 - (A0_inv * FX0);
    vector<double> X1_vec = vector_from_matrix(X1);

    // step 5
    vector<AD> FX1_AD_vec;         
    vector<double> FX1_vec;       
    matrix FX1(no_of_elements, 1); 

    while (error > epsilon)
    {
        // step 6
        FX1_AD_vec = gen_PDE.evaluate_PDE(X1_vec);
        FX1_vec = makeDoubleVector(FX1_AD_vec);
        FX1 = matrix_from_vector(FX1_vec, no_of_elements, 1);

        matrix y1 = FX1 - FX0;
        matrix s1 = X1 - X0;

        // step 7
        matrix s1t = s1.transpose();
        matrix temp1 = s1t * A0_inv; 
        matrix temp2 = temp1 * y1;
        double s1t_x_A0_inv_x_y1 = temp2.M[0][0]; 

        // step 8
        matrix A1_inv = A0_inv + ((((s1 - (A0_inv * y1)) * s1t) * A0_inv) / (s1t_x_A0_inv_x_y1));

        // step 9
        matrix X2 = X1 - (A1_inv * FX1);

        // step 10
        vector<double> X2_vec = vector_from_matrix(X2);
        double error = norm(X2_vec, X1_vec);

        for (int i = 0; i < no_of_elements; i++)
        {
            for (int j = 0; j < no_of_elements; j++)
            {
                A0_inv.M[i][j] = A1_inv.M[i][j];
            }
            X0.M[i][0] = X1.M[i][0];
            X1.M[i][0] = X2.M[i][0];
            FX0.M[i][0] = FX1.M[i][0];
        }

        cout << "X2 = " << endl;
        X2.print();
        cout << endl;

        X1_vec = X2_vec;
        count++;
    }
    vector<double> final_soln = X1_vec;

    for (int k = 0; k < gen_PDE.domain.ind_var_Points.size(); k++)
    {
        AD temp(final_soln[k], k, gen_PDE.domain.ind_var_Points.size());
        gen_PDE.domain.ind_var_Points[k].set_u(temp);
    }
    return gen_PDE.domain;
}



#endif