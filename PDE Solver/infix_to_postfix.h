#ifndef INFIX_TO_POSTFIX
#define INFIX_TO_POSTFIX
#include <iostream>
using namespace std;
#include <cmath>
#include <string>
#include <vector>
#include "AD.h"
#define MAX 1000

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

#endif