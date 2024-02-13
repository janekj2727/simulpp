// g++ -Wall test_rand.cpp Vector.cpp -o test_rand -lm -I ../incl

/*
*  Main for testing random number generators in math_utils
*
*  Author: JJ
*  Date:
*
*/

#include "math_utils.hpp"
#include "Vector.hpp"
#include <iostream>
#include <cstdio>

int main(int argc, char **argv)
{
    // std::cout << "uniform: " << randuni() << "\n";
    // std::cout << "gaussian: " << randgauss() << "\n";

    // int i;
    // // double u1, u2;
    // for (i = 0; i < 1000000; i++)
    // {
    //     std::cout << randuni() << "\n";
    //     // u1 = (double) generator() / generator.max();
    //     // u2 = (double) generator() / generator.max();
    //     // std::cout << u1 << "    " << sqrt(-2*log(u1))*cos(2*M_PI*u2) << "\n";
    // }

    Vector a(stdin);
    double stddev;


    std::cout << a << "\n" << a.Mean() << "\n" << a.MeanWithError(stddev, false) << " " << stddev << "\n"; 
    // std::cout << CalculateAngle(a, b) << std::endl;
    return 0;
}

