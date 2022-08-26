#ifndef DISCRETIZER_H
#define DISCRETIZER_H

#include <iostream>
#include <fstream>
#include "AD.h"
#include <vector>
#include <cmath>
#include <string>
using namespace std;
#define MAX 1000

//______________________________________    Infix to Postfix    ________________________________________//

class stringList
{
private:
    string elements[MAX];

public:
    stringList();
    stringList(int);

    stringList infix_to_postfix(string);
    AD evaluate(double, double, int, AD);

    void display();
    // void pushElement(string);
    void setElement(int, string);
    string getElement(int);

    bool is_var_char(char c);
    bool is_num_char(char c);
    bool isOperator(string s);
    int precedence(string s);
    bool associativity(string s);
    bool isMathFunction(string fname);
    vector<string> tokenize(string);
};

template <class T>
class Stack
{
private:
    T elements[MAX];
    int size;

public:
    Stack();
    bool isEmpty();
    bool isFull();
    int getTop();
    T getTopElement();
    void push(T);
    T pop();
};

// ________________________________________  Class StringList Definitions  _______________________________________//

stringList::stringList()
{
    for (int i = 0; i < MAX; i++)
    {
        elements[i] = "";
    }
}

void stringList::display()
{
    for (int i = 0; elements[i] != ""; i++)
    {
        if (elements[i] != "$")
        {
            cout << elements[i];
        }
    }
    cout << endl;
}

void stringList::setElement(int i, string s)
{
    this->elements[i] = s;
}

string stringList::getElement(int i)
{
    return elements[i];
}

//__________________________<<<<< stringList : infix_to_postfix Function Definitons >>>>>________________________//

stringList stringList::infix_to_postfix(string expression)
{
    vector<string> postfix_vec(0);
    Stack<string> sStack;
    int n = expression.length();
    vector<string> infix_vec = tokenize(expression);

    for (int i = 0; i < infix_vec.size(); i++)
    {
        string s = infix_vec[i];
        // If an operand is scanned then add it to the postfixExpression
        if ((is_var_char(s[0]) && !isMathFunction(s)) || (is_num_char(s[0])) || s == "$")
        {
            postfix_vec.push_back(s);
        }
        // If the scanned item is '(' push it into the sStack
        else if (s == "(" || isMathFunction(s) || s == "#")
        {
            sStack.push(s);
        }
        // If the scanned item is ')' pop the Stack and append it to the PFE until top element becomes '(' then pop ')'.
        else if (s == ")")
        {
            while (sStack.getTopElement() != "(" && !sStack.isEmpty())
            {
                postfix_vec.push_back(sStack.getTopElement());
                string sStackTop = sStack.pop();
            }
            string popped_openParanthesis = sStack.pop();
            if (sStack.getTopElement() == "#" && !sStack.isEmpty())
            {
                postfix_vec.push_back(sStack.getTopElement()); // add the sStack top element to the postfixExpression
                string sStackTop_hash = sStack.pop();
                postfix_vec.push_back(sStack.getTopElement()); // add the sStack top element to the postfixExpression
                string sStackTop_math_func_name = sStack.pop();
            }
        }
        // If an operator is scanned
        else if (isOperator(s))
        {
            if (sStack.isEmpty())
            {
                sStack.push(s);
            }
            else if (precedence(s) > precedence(sStack.getTopElement()))
            { // push s into sStack;
                sStack.push(s);
            }
            else if (precedence(s) == precedence(sStack.getTopElement()))
            { // check associativity
                if (associativity(s))
                { // if ltr associativity, then first add and pop the sStack top element to PFE and then push "s" into sStack.
                    postfix_vec.push_back(sStack.getTopElement());
                    string op = sStack.pop();
                    sStack.push(s);
                }
                // if associativity is rtl. in case of '^' operator.
                else
                { // simply push the 's' into the sStack.
                    sStack.push(s);
                }
            }
            // if the precedence(s)  is less than or equal precedence(sStack Top element)
            else
            {
                while (!sStack.isEmpty() && (precedence(s) <= precedence(sStack.getTopElement())))
                { // while it holds true add the top element of sStack to PFE and then pop it from sStack.
                    postfix_vec.push_back(sStack.getTopElement());
                    string pop_s = sStack.pop();
                }
                sStack.push(s); // finally push 's' here in this case the lower precedence operator
            }
        }
    }
    // while the sStack isn't empty, pop and add the elements of it to the PFE.
    while (!sStack.isEmpty())
    {
        postfix_vec.push_back(sStack.getTopElement());
        string pop_d = sStack.pop();
    }

    stringList postfix;
    for (int i = 0; i < postfix_vec.size(); i++)
    {
        postfix.setElement(i, postfix_vec[i]);
    }

    return postfix;
}

//________________________________<<<<< stringList : Evaluate Function Definitons >>>>__________________________//

AD stringList ::evaluate(double xValue, double yValue, int no_of_ind_vars, AD U)
{
    // x$y$*x$u$12.3$/+#sin+
    Stack<AD> adStack;
    for (int i = 0; elements[i] != ""; i++)
    {
        string s = elements[i];

        if (is_num_char(s[0]))
        {
            double num = stod(s);
            AD constant_number(num, no_of_ind_vars);
            adStack.push(constant_number);
        }
        else if (s == "x" || s == "y")
        {
            if (s == "x")
            {
                AD x(xValue, no_of_ind_vars); // x will behave as constant in context of discrete points.
                adStack.push(x);
            }
            else if (s == "y")
            {
                AD y(yValue, no_of_ind_vars); // y will behave as constant in context of discrete points.
                adStack.push(y);
            }
        }
        else if (s == "u")
        {
            adStack.push(U);
        }
        else if (isOperator(s))
        {
            AD var2, var1, res;
            if (s == "+")
            {
                var2 = adStack.getTopElement();
                adStack.pop();
                var1 = adStack.getTopElement();
                adStack.pop();
                res = var1 + var2;
            }
            else if (s == "-")
            {
                var2 = adStack.getTopElement();
                adStack.pop();
                if(adStack.isEmpty())
                {   AD minus_ONE(-1,no_of_ind_vars);
                    res = minus_ONE * var2;
                }   
                else {
                    var1 = adStack.getTopElement();
                    adStack.pop();
                    res = var1 - var2;
                }
            }
            else if (s == "*")
            {
                var2 = adStack.getTopElement();
                adStack.pop();
                var1 = adStack.getTopElement();
                adStack.pop();
                res = var1 * var2;
            }
            else if (s == "/")
            {
                var2 = adStack.getTopElement();
                adStack.pop();
                var1 = adStack.getTopElement();
                adStack.pop();
                res = var1 / var2;
            }
            else if (s == "^")
            {
                var2 = adStack.getTopElement();
                adStack.pop();
                var1 = adStack.getTopElement();
                adStack.pop();
                res = var1 ^ var2;
            }
            adStack.push(res);
        }
        else if (isMathFunction(s))
        { // x$y$*x$u$12.3$/+#sin+
            AD var = adStack.getTopElement();
            adStack.pop();
            AD res;
            if (s == "sin")
            {
                res = sin(var);
            }
            else if (s == "cos")
            {
                res = cos(var);
            }
            else if (s == "tan")
            {
                res = tan(var);
            }
            else if (s == "cosec")
            {
                res = cosec(var);
            }
            else if (s == "sec")
            {
                res = sec(var);
            }
            else if (s == "cot")
            {
                res = cot(var);
            }
            else if (s == "arcsin")
            {
                res = arcsin(var);
            }
            else if (s == "arccos")
            {
                res = arccos(var);
            }
            else if (s == "arctan")
            {
                res = arctan(var);
            }
            else if (s == "sinh")
            {
                res = sinh(var);
            }
            else if (s == "cosh")
            {
                res = cosh(var);
            }
            else if (s == "tanh")
            {
                res = tanh(var);
            }
            else if (s == "log")
            {
                res = log(var);
            }
            else if (s == "exp")
            {
                res = exp(var);
            }
            else if (s == "abs")
            {
                res = abs(var);
            }
            adStack.push(res);
        }
    }
    AD answer = adStack.getTopElement();
    return answer;
}

//________________________________<<<<< stringList : Friend Function Definitons >>>>>___________________________//

bool stringList::is_var_char(char c)
{
    if (('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z') || c == '_')
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool stringList ::is_num_char(char c)
{
    if (('0' <= c && c <= '9') || c == '.')
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool stringList::isOperator(string s)
{ // check if "s" is an operator or not.
    if (s == "^" || s == "*" || s == "/" || s == "+" || s == "-")
        return true;
    else
        return false;
}

int stringList ::precedence(string s)
{ // assign precedence to operators.
    if (s == "^")
        return 3;
    else if (s == "*" || s == "/")
        return 2;
    else if (s == "+" || s == "-")
        return 1;
    else
        return 0;
}

bool stringList::associativity(string s)
{
    // returns "true" for left to right associativity and "false" for right to left associativity.
    if (s == "*" || s == "/" || s == "+" || s == "-") // ltr associativity
        return true;
    else if (s == "^") // rtl associativity
        return false;
}

bool stringList ::isMathFunction(string fname)
{
    // checks if 'c' is a MathFunction. S for sin, C for cos and so on, look into
    // replace_Math_Functions_Name_By_Char fuction for more information on char codes.
    if (fname == "sin" || fname == "cos" || fname == "tan" ||
        fname == "cosec" || fname == "sec" || fname == "cot" ||
        fname == "arcsin" || fname == "arccos" || fname == "arctan" ||
        fname == "sinh" || fname == "cosh" || fname == "tanh" ||
        fname == "log" || fname == "exp" || fname == "abs")
    {
        return true;
    }
    else
    {
        return false;
    }
}

vector<string> stringList ::tokenize(string input)
{
    int inp_len = input.length();

    vector<string> infix_vector;

    for (int i = 0; i < inp_len; i++)
    {
        string opr(1, input[i]);
        if (isOperator(opr) || input[i] == '(' || input[i] == ')')
        {
            string temp(1, input[i]);
            infix_vector.push_back(temp);
        }
        else if (is_num_char(input[i]))
        {
            string number = "";
            while (is_num_char(input[i]) && i < inp_len)
            {
                number = number + input[i];
                // cout << number << endl;
                i++;
            }
            i--;
            infix_vector.push_back(number);
            infix_vector.push_back("$");
        }
        else if (is_var_char(input[i]))
        {

            string word = "";
            while ((is_var_char(input[i]) || is_num_char(input[i])) && input[i] != '.' && i < inp_len)
            {
                word = word + input[i];
                i++;
            }
            i--;
            if (isMathFunction(word))
            {
                infix_vector.push_back(word);
                infix_vector.push_back("#");
            }
            else
            {
                infix_vector.push_back(word);
                infix_vector.push_back("$");
            }
        }
    }

    return infix_vector;
}

// --------------------------------------------------------------------------------------------------------------//

// ________________________________________  Template Stack Definitions  ________________________________________//
template <class T>
Stack<T>::Stack()
{
    size = 0;
}

template <class T>
bool Stack<T>::isEmpty()
{
    if (size == 0)
    {
        return true;
    }
    return false;
}

template <class T>
bool Stack<T>::isFull()
{
    if (size >= MAX)
    {
        return true;
    }
    return false;
}

template <class T>
int Stack<T>::getTop()
{
    return size;
}

template <class T>
T Stack<T>::getTopElement()
{
    return elements[size - 1];
}

template <class T>
void Stack<T>::push(T value)
{
    size++;
    elements[size - 1] = value;
}

template <class T>
T Stack<T>::pop()
{
    size--;
    return elements[size];
}

//-------------------------------------------------------------------------------------------------------//
class Discrete_Point
{
private:
    int i, j;
    double x_val, y_val;
    AD u;

public:
    Discrete_Point();
    Discrete_Point(int, int, double, double);
    void set_i(int);
    void set_j(int);
    void set_u(AD);
    void set_x_val(double);
    void set_y_val(double);
    int get_i();
    int get_j();
    AD get_u();
    double get_x_val();
    double get_y_val();
};



//____________________________________   Class Discrete_Point Definitions   _________________________________________//

Discrete_Point::Discrete_Point()
{
    this->i = 0;
    this->j = 0;
    this->x_val = 0;
    this->y_val = 0;
}

Discrete_Point ::Discrete_Point(int i, int j, double x_val, double y_val)
{
    this->i = i;
    this->j = j;
    this->x_val = x_val;
    this->y_val = y_val;
    AD u();
}

void Discrete_Point ::set_i(int i)
{
    this->i = i;
}

void Discrete_Point ::set_j(int j)
{
    this->j = j;
}

void Discrete_Point ::set_u(AD u)
{
    this->u = u;
}

void Discrete_Point ::set_x_val(double x_val)
{
    this->x_val = x_val;
}

void Discrete_Point ::set_y_val(double y_val)
{
    this->y_val = y_val;
}

int Discrete_Point ::get_i()
{
    return this->i;
}

int Discrete_Point ::get_j()
{
    return this->j;
}

AD Discrete_Point ::get_u()
{
    return this->u;
}

double Discrete_Point ::get_x_val()
{
    return this->x_val;
}

double Discrete_Point ::get_y_val()
{
    return this->y_val;
}

//---------------------------------------------------------------------------------------------------------------//

//____________________________________   Class Discretize Definitions   _________________________________________//

class Discretize
{
private:
    double a, b, c, d;

public:
    int m, n;
    double x_value, y_value;
    double delta_x, delta_y;
    vector<Discrete_Point> ind_var_Points;
    vector<Discrete_Point> baoundary_Points;
    vector<double> F;
    Discretize();
    Discrete_Point get_Point(int, int);
    Discretize(double, double, double, double, int, int, string, string, string, string);
    void display();
};

Discretize ::Discretize()
{
    this->a = 0;
    this->b = 0;
    this->c = 0;
    this->d = 0;
    this->m = 0;
    this->n = 0;
}

Discretize ::Discretize(double a, double b, double c, double d, int m, int n, string f1_y, string f2_y, string g1_x, string g2_x)
{
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
    this->m = m;
    this->n = n;
    this->F = vector<double>((m - 1) * (n - 1));

    delta_x = abs((b - a) / m);
    delta_y = abs((d - c) / n);

    stringList infix_f1_y, infix_f2_y, infix_g1_x, infix_g2_x;
    stringList postfix_f1_y, postfix_f2_y, postfix_g1_x, postfix_g2_x;
    postfix_f1_y = infix_f1_y.infix_to_postfix(f1_y);
    postfix_f2_y = infix_f2_y.infix_to_postfix(f2_y);
    postfix_g1_x = infix_g1_x.infix_to_postfix(g1_x);
    postfix_g2_x = infix_g2_x.infix_to_postfix(g2_x);

    for (int i = 1; i < m; i++)
    {
        for (int j = 1; j < n; j++)
        {
            Discrete_Point interior_point(i, j, a + delta_x * i, c + delta_y * j);
            double initial_guess = 0.1;
            AD u(initial_guess, (n - 1) * (i - 1) + j - 1, (m - 1) * (n - 1));
            interior_point.set_u(u);
            ind_var_Points.push_back(interior_point);
        }
    }
    AD zero(0, (m - 1) * (n - 1));
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            x_value = a + delta_x * i;
            y_value = c + delta_y * j;

            if (i == 0 || (i == 0 && j == 0))
            {
                Discrete_Point bPoint(i, j, x_value, y_value);
                AD f1_y_val = postfix_f1_y.evaluate(x_value, y_value, (m - 1) * (n - 1), zero);
                bPoint.set_u(f1_y_val);
                baoundary_Points.push_back(bPoint);
            }
            else if (i == m || (i == m && j == m))
            {
                Discrete_Point bPoint(i, j, x_value, y_value);
                AD f2_y_val = postfix_f2_y.evaluate(x_value, y_value, (m - 1) * (n - 1), zero);
                bPoint.set_u(f2_y_val);
                baoundary_Points.push_back(bPoint);
            }
            else if (j == 0 && i != 0)
            {
                Discrete_Point bPoint(i, j, x_value, y_value);
                AD g1_x_val = postfix_g1_x.evaluate(x_value, y_value, (m - 1) * (n - 1), zero);
                bPoint.set_u(g1_x_val);
                baoundary_Points.push_back(bPoint);
            }
            else if (j == n && i != n)
            {
                Discrete_Point bPoint(i, j, x_value, y_value);
                AD g2_x_val = postfix_g2_x.evaluate(x_value, y_value, (m - 1) * (n - 1), zero);
                bPoint.set_u(g2_x_val);
                baoundary_Points.push_back(bPoint);
            }
            else
            {
                Discrete_Point bPoint(i, j, x_value, y_value);
                bPoint.set_u(zero);
                baoundary_Points.push_back(bPoint);
            }
        }
    }
}

Discrete_Point Discretize ::get_Point(int i, int j)
{
    if (i == 0 || i == m || j == 0 || j == n)
    {
        return this->baoundary_Points[(n + 1) * i + j];
    }
    else
    {
        return this->ind_var_Points[(n - 1) * (i - 1) + (j - 1)];
    }
}

void Discretize::display()
{
    
    cout << "X value" << setw(15) << "Y value" << setw(15) << "Z(X,Y) value"<< endl; 
    for(int i = 0; i < m; i++)
    {
        for(int j =0; j < n; j++)
        {
        cout.precision(5);
        cout << get_Point(i,j).get_x_val() << setw(15) << get_Point(i,j).get_y_val() << setw(15) << get_Point(i,j).get_u().getf()<< endl; 
        }
    }
}
//---------------------------------------------------------------------------------------------------------------//
void create_csv_file(Discretize Solution)
{   //
    ofstream myFile;
    myFile.open("plot_data.csv");

    for (int i = 0; i <=Solution.m; i++)
    {
        for (int j = 0; j <=Solution.n; j++)
        {
            myFile << Solution.get_Point(i, j).get_x_val() << "," << Solution.get_Point(i, j).get_y_val() << "," << Solution.get_Point(i, j).get_u().getf() << endl;
        }
    }
    myFile.close();
}

#endif