/*
 * definitions of small math utilities – factorial, greatest common divisor,...
 * originally for simul++
 * Author: JJ
 * Date: Apr 2021
 */

#include "math_utils.hpp"
#include <chrono>
#include <cmath>
#include <random>
#include <cassert>
#include "Matrix.hpp"

// #include <iostream>

// factorial for predictor matrix
int fact(int x)
{
    int y = 1;
    int i;
    for (i = 1; i < x + 1; i++)
    {
        y *= i;
    }
    return y;
}

// calculate greatest common divisor of two integers a and b
int gcd(int a, int b)
{
    if (a == 0 || b == 0)
        return 0;
    else if (a == b)
        return a;
    else if (a > b)
        return gcd(a - b, b);
    else
        return gcd(a, b - a);
}

// calculate greatest common divisor of three integers in vector abc (zeros ignored)
// return 0 only if all 3 ints == 0
int gcd3(int *abc)
{
    int orig[3];
    orig[0] = abc[0];
    orig[1] = abc[1];
    orig[2] = abc[2];
    int d, e;
    if (abc[0] == 0)
    {
        abc[0] = (abc[1] > abc[2]) ? abc[1] : abc[2];
    }
    if (abc[1] == 0)
    {
        abc[1] = (abc[0] > abc[2]) ? abc[0] : abc[2];
    }
    if (abc[2] == 0)
    {
        abc[2] = (abc[0] > abc[1]) ? abc[0] : abc[1];
    }
    d = gcd(abc[0], abc[1]);
    e = abc[2];
    abc[0] = orig[0];
    abc[1] = orig[1];
    abc[2] = orig[2];
    return gcd(e, d);
}

// calculate greatest common divisor of three integers in vector abc (zeros ignored)
// return 0 only if all 3 ints == 0
int gcd4(int *abcd)
{
    int orig[4];
    orig[0] = abcd[0];
    orig[1] = abcd[1];
    orig[2] = abcd[2];
    orig[3] = abcd[3];
    int e, f;
    if (abcd[0] == 0)
    {
        abcd[0] = (abcd[1] > abcd[2]) ? abcd[1] : abcd[2];
    }
    if (abcd[1] == 0)
    {
        abcd[1] = (abcd[0] > abcd[2]) ? abcd[0] : abcd[2];
    }
    if (abcd[2] == 0)
    {
        abcd[2] = (abcd[0] > abcd[1]) ? abcd[0] : abcd[1];
    }
    if (abcd[3] == 0)
    {
        abcd[3] = (abcd[0] > abcd[1]) ? abcd[0] : abcd[1];
    }
    f = gcd(abcd[0], abcd[1]);
    e = gcd(abcd[2], abcd[3]);
    abcd[0] = orig[0];
    abcd[1] = orig[1];
    abcd[2] = orig[2];
    return gcd(e, f);
}

// random numbers equally distributed (0,1)
double randuni()
{
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::ranlux48 generator(seed);
    static uint_fast64_t max = generator.max();

    // std::cout << "min: " << generator.min() << "\n";

    return (double)generator() / max;
}

// random numbers with gaussian distribution (µ = 0, sigma = 1.0)
double randgauss()
{
    static std::normal_distribution<double> normal(0.0, 1.0);
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::ranlux48 generator(seed);

    return normal(generator);
}

// all-pairs-shortest-path in graph given by connectivity matrix A
// all edges are bidirectional and unweighted
// primarily for bond-distance in a Molecule
// Seidel algorithm accorging to J. Comp. Sys. Sci. 51, 400–403 (1995) (doi:10.1006/jcss.1995.1078)
// connectivity matrix must be symmetrized, otherwise returns a mess
Matrix *shortestpath(const Matrix A, int recursionlevel)
{
    Matrix *Z, *B, *D, *T, *X;
    int Nrow, Ncol;
    int i, j;
    bool done = false;

    Nrow = A.GetNumberOfRows();
    Ncol = A.GetNumberOfCols();
    assert(Nrow == Ncol);
    Z = new Matrix(Nrow, Ncol);
    B = new Matrix(Nrow, Ncol);
    D = new Matrix(Nrow, Ncol);
    X = new Matrix(Nrow, Ncol);

    (*Z) = A * A;
    for (i = 0; i < Nrow; i++)
    {
        for (j = 0; j < Ncol; j++)
        {
            if (i == j)
            {
                B->operator()(i, j) = 0.0;
            }
            else if ((A.Read(i, j) == 1.0) || (Z->operator()(i, j) > 0.0))
            {
                B->operator()(i, j) = 1.0;
            }
        }
    }
    done = true;
    for (i = 0; i < Nrow; i++)
    {
        for (j = 0; j < Ncol; j++)
        {
            if ((B->operator()(i, j) == 0.0) && (i != j))
            {
                done = false;
            }
        }
    }
    if (done || (recursionlevel > Nrow)) // prevents infinity recursion in case of wrong input
    {
        (*D) = (*B) * 2.0 - A;
        delete Z;
        delete B;
        delete X;
        return D;
    }
    else
    {
        T = shortestpath(*B, recursionlevel + 1); // safety check
        (*X) = (*T) * A;
        for (i = 0; i < Nrow; i++)
        {
            for (j = 0; j < Ncol; j++)
            {
                if (X->Read(i, j) >= T->Read(i, j) * Z->Read(j, j))
                {
                    D->operator()(i, j) = 2.0 * T->Read(i, j);
                }
                else
                {
                    D->operator()(i, j) = 2.0 * T->Read(i, j) - 1.0;
                }
            }
        }
        delete Z;
        delete B;
        delete X;
        delete T;
        return D;
    }

    // You cannot reach this
    delete Z;
    delete B;
    delete X;
    return D;
}