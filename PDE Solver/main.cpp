#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include "AD.h"
#include "Matrix.h"
// #include "infix_to_postfix.h"
#include "Discretizer.h"
#include "NonLinearSolver.h"
#include "LinearSolver.h"
int main()
{
    // // General form of the nonlinear PDE :    A(x,y,u)*d2u/dx2 + B(x,y,u)*d2u/dy2 + C(x,y,u)*du/dx + D(x,y,u)*du/dy + E(x,y,u) = 0 ;
    // // on the range (x * y) where x -> (a,b) and y -> (c,d) ;
    // // x,y -> independent variable,    u -> dependent variable ;

    // //_________________________________________________ Take the input from the user __________________________________________________//

    // ask the user for coeffs of the nonlinear PDE
    cout << endl;
    cout << "General form of the nonlinear PDE : " << endl
         << endl;
    cout << "A(x,y,u)*d2u/dx2 + B(x,y,u)*d2u/dy2 + C(x,y,u)*du/dx + D(x,y,u)*du/dy + E(x,y,u) = 0 " << endl
         << endl;

    cout << "Set the coefficients of the nonlinear PDE" << endl;
    string A, B, C, D, E;
    cout << "Enter A(x,y,u) : ";
    cin >> A;
    cout << "Enter B(x,y,u) : ";
    cin >> B;
    cout << "Enter C(x,y,u) : ";
    cin >> C;
    cout << "Enter D(x,y,u) : ";
    cin >> D;
    cout << "Enter E(x,y,u) : ";
    cin >> E;
    cout << endl;

    // ask the user for the domain [a,b] x [c,d] values
    cout << "Set domain values" << endl;
    cout << "Range of x = [a,b] and Range of y = [c,d] , where a, b, c, d are natural numbers" << endl;
    double a, b, c, d;
    cout << "Enter a : ";
    cin >> a;
    cout << "Enter b : ";
    cin >> b;
    cout << "Enter c : ";
    cin >> c;
    cout << "Enter d : ";
    cin >> d;
    cout << endl;
    if ((a > b) && (c > d))
    {
        cout << "Please make sure that b > a and d > c " << endl;
    }

    // ask the user for boundary conditions
    cout << "Set boundary conditions" << endl;
    cout << "u(a,y) = f1(y),    u(b,y) = f2(y),    u(x,c) = g1(x),    u(x,d) = g2(x)" << endl;
    string f1_y, f2_y, g1_x, g2_x;
    cout << "Enter f1(y) : ";
    cin >> f1_y;
    cout << "Enter f2(y) : ";
    cin >> f2_y;
    cout << "Enter g1(x) : ";
    cin >> g1_x;
    cout << "Enter g2(x) : ";
    cin >> g2_x;
    cout << endl;

    // ask the user for the divisions along x axis (m) and divisions along y axis (n);
    int m, n;
    cout << "Enter the number of divisions along X axis : ";
    cin >> m;
    cout << "Enter the number of divisions along Y axis : ";
    cin >> n;
    cout << endl;

    // ask the user for to choose from the available solver types;
    int NonLinearChoice, LinearChoice;
    cout << "Choose Non Linear Solver" << endl; 
    cout << "[1] for Newton Solver" << endl; 
    cout << "[2] for Broyden Solver" << endl;
    cin >> NonLinearChoice;
    cout<< endl;
    cout << "Choose Linear Solver" << endl; 
    cout << "[1] for Gauss Elimination Solver" << endl; 
    cout << "[2] for LU Decomposition Solver" << endl; 
    cout << "[3] for TriDiagonal Solver" << endl; 
    cout << "[4] for Gauss Jacobi Solver" << endl; 
    cout << "[5] for Gauss Seidal Solver" << endl; 
    cout << "[6] for Successive Over Relaxation Solver" << endl; 
    cin >> LinearChoice;
    cout<< endl;

    // //___________________________________________    Perform the Operations on input     ____________________________________________//

    Discretize domain(a, b, c, d, m, n, f1_y, f2_y, g1_x, g2_x);
    PDE gen_PDE(A, B, C, D, E, domain);
    Discretize ans = Solve(gen_PDE,NonLinearChoice,LinearChoice);
    ans.display();
    // create_csv_file(ans);        // Creates a CSV file containing (X,Y,Z) data. 

    return 0;
}

