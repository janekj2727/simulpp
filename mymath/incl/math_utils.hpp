/*
* definitions of small math utilities – factorial, greatest common divisor,...
* originally for simul++
* Author: JJ
* Date: Apr 2021
*/

#ifndef MATHUTILSHEADER
#define MATHUTILSHEADER

class Matrix;

// factorial for predictor matrix
int fact(int x);

// calculate greatest common divisor of two integers a and b
int gcd(int a, int b);

// calculate greatest common divisor of three integers in vector abc (zeros ignored)
// return 0 only if all 3 ints == 0
int gcd3(int *abc);

// calculate greatest common divisor of four integers in vector abc (zeros ignored)
int gcd4(int *abc);

// random numbers equally distributed (0,1) (uniform)
double randuni();

// random numbers with gaussian distribution (µ = 0, sigma = 1.0)
double randgauss();

// all-pairs-shortest-path solution for graph given by matrix A
Matrix *shortestpath(Matrix A, int recursionlevel = 0);

#endif