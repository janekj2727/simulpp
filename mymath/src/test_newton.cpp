/* Test of Newton Method*/

#include <cmath>
#include <iostream>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "NewtonMethod.hpp"
#include "general_utils.hpp"
#include "MDTable.hpp"

Vector testFunction(Vector& x);

int main(int argc, char **argv)
{
    NewtonMethod *method;
    Vector x0(2), res(2);
    Vector (*function)(Vector&);

    function = &testFunction;

    method = new NewtonMethod(function, x0, 30);
    res = method->Solve();

    std::cout << "Result is:\n" << res << "\n";

    // Matrix *A;

    // A = new Matrix(2,2);

    // A->operator()(0,0) = 2.0;
    // A->operator()(1,0) = 3.0;
    // A->operator()(0,1) = 1.0;
    // std::cout << "Matrix before: \n" << *A << "\n";
    // A->operator()(1,0) *= 2.0;
    // std::cout << "Matrix after A(1,0) *= 2.0: \n" << *A << "\n";

    // std::cout << "Vector before: " << b << "\n";
    // b(1) -= 2.0;
    // std::cout << "Vector after b(0) -= 2.0: " << b << "\n";

    // MDTable<std::string, double, int, std::string> vt({"Name", "Weight", "Age", "Brother"},
    //                                                         10);

    // vt.setColumnFormat({MDTableColumnFormat::AUTO, MDTableColumnFormat::SCIENTIFIC, MDTableColumnFormat::AUTO, MDTableColumnFormat::AUTO});
    // vt.setColumnPrecision({0, 12, 3, 0});

    // vt.addRow("Cody", 180.26885448874, 40, "John");
    // vt.addRow("David", 1750.3, 38, "Andrew");
    // vt.addRow("Robert", 14000.3, 27, "Fande");

    // vt.print(std::cout);

    // Matrix M(2, 2);
    // M(0,0) = 0.0;
    // M(0,1) = -23.436536454515655;
    // M(1, 0) = -1.112548564458774545e-10;
    // M(1, 1) = 18.5;

    // std::vector<std::string> colh({"Prvni nadpis", "Druhy dlooouhy nadpis"});
    // std::vector<std::string> rowh({"Prvni radek", "Druhy"});

    // std::cout << "Matrix M by operator <<:\n" << M << "\n";
    // std::cout << "Matrix M by PrintMDTable:\n" << std::setprecision(8);
    // M.PrintMDTable(std::cout, colh, rowh);


    return 0;
}

Vector testFunction(Vector& x)
{
    Vector result(2);

    if (x.GetSize() != 2)
    {
        std::cerr << "ERROR: Vector passed to 'testFunction' has wrong size!!!\n";
        return result;
    }

    result[0] = x[0] + x[0] * x[1] - 4;
    result[1] = x[0] + x[1] - 3;

    return result;
}
//no instance of constructor "NewtonMethod::NewtonMethod" matches the argument list -- argument types are: (Vector (*)(Vector &x), Vector)