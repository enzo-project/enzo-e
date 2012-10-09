// File:    test.cpp
// Author:  Brian Allen Vanderburg II
// Purpose: ExprEval 3.x test program
//------------------------------------------------------------------------------

// Includes
#include <iostream>
#include <ctime>

#include "../expreval.h"

using namespace std;
using namespace ExprEval;

    
class TestNode : public FunctionNode
{
public:
    TestNode(Expression *expr) : FunctionNode(expr)
    {
        SetArgumentCount(0, 0, 0, 0);
    }
        
    double DoEvaluate()
    {
        return 123.456;    
    }
};
    
class TestFactory : public FunctionFactory
{
public:
    string GetName() const
    {
        return "test";
    }
        
    FunctionNode *DoCreate(Expression *expr)
    {
        return new TestNode(expr);
    }
};        
    
void HandleException(Exception &e)
{
    switch(e.GetType())
    {
        case Exception::Type_SyntaxException:
            cout << "Syntax Error\n";
            break;
            
        case Exception::Type_EmptyExpressionException:
            cout << "Empty Expression\n";
            break;
            
        default:
            cout << "Other Error\n";
            break;        
    }
}

int main()
{
    string expr;
    long pos, count;
    
    cout << "Enter expression: ";
    cin >> expr;
    
    
    cout << "Enter count: ";
    cin >> count;
    
    try
    {
        ValueList v;
        FunctionList f;
        Expression e;
        
        v.AddDefaultValues();
        f.AddDefaultFunctions();
        f.Add(new TestFactory());

        e.SetValueList(&v);
        e.SetFunctionList(&f);

        e.Parse(expr);
        
        double result;
        
        // Start a timer
        time_t start = time(NULL);
        
        for(pos = 0; pos < count; pos++)
        {
            result = e.Evaluate();
        }
            
        // Determine total time
        double total = difftime(time(NULL), start);
        
        cout << "Total time: " << total << endl;
        if(total != 0.0)
        {
            cout << "Average/sec: " << (double)count / total << endl;
        }
            
        // Variable dump
        cout << "Value dump" << endl;
        cout << "-------------------------" << endl;
        
        ValueList::size_type i, c;
        c = v.Count();
        for(i = 0; i < c; i++)
        {
            string name;
            double value;
            
            v.Item(i, &name, &value);
            cout << name.c_str() << " = " << value << endl;
        }
            
        cout << endl;
    }
    catch(Exception &e)
    {
        HandleException(e);
    }
}
