// This file was created automatically with /home/janekj/Programs/simul++/bin/mymath.sh
#include "math_utils.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "ARMA.hpp"
#include "NewtonMethod.hpp"
#include "LinearSystem.hpp"
#include "LinearInterpolator.hpp"
#include "Linear3PInterpolator.hpp"
#include "Linear4PInterpolator.hpp"
#include "HermiteCubSplines.hpp"
#include "NaturalCubSplines.hpp"
#include "MacsimusQuadSplines.hpp"
#include "MacsimusHyperbSplines.hpp"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <ctime>
#include <chrono>

int main(int argc, char **argv)
{
std::vector<double> knots;
std::vector<double> derivatives, derivatives2, derivatives3;
int i;
double h = 1.0/2.0;
int N = 30;
double max;
double r2, r;
double alpha = 0.4;
double min = 1.0;

for (i=0; i<(N+1); i++)
{
r2 = min + i * h;
r = sqrt(r2);
knots.push_back(erfc(alpha * r)/r);
derivatives.push_back(((-2.0 * alpha * exp(-alpha*alpha*r2))/(sqrt(M_PI)*r) - erfc(alpha*r)/r2)/(2 * r));
}

derivatives2.push_back(*(derivatives.begin()));
derivatives2.push_back(*(derivatives.end() - 1));
derivatives3.push_back(*(derivatives.end() - 1));


HermiteCubSplines<double> hermite(knots, min + N*h, min, derivatives);
NaturalCubSplines<double> natural(knots, min + N*h, min, derivatives2);
NaturalCubSplines<double> natural2(knots, min + N*h, min);
MacsimusQuadSplines<double> quadratic(knots, min + N*h, min, derivatives3);
MacsimusHyperbSplines<double> hyperbolic(knots, min + N*h, min, derivatives3);
LinearInterpolator<double> linear(knots, min + N*h, min);
Linear3PInterpolator<double> linear3(knots, min + N*h, min);
Linear4PInterpolator<double> linear4(knots, min + N*h, min);

if (hermite.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (natural.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (natural2.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (quadratic.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (hyperbolic.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (linear.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (linear3.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";
if (linear4.ExtendToZero() < 0) std::cerr << "Not able to extend to zero!\n";

std::vector<double> xval;
std::vector<double> xsqrt;
std::vector<double> yher;
std::vector<double> ynat;
std::vector<double> ynat2;
std::vector<double> yquad;
std::vector<double> yhyp;
std::vector<double> ylin;
std::vector<double> ylin3;
std::vector<double> ylin4;
std::vector<double> yexact;

max = (N * h);
h = 1.0/1024.0;

for (i=1025; i<floor(max/h) - 1025; i++)
{
r2 = min + i * h;
r = sqrt(r2);
yexact.push_back((((-2.0 * alpha * exp(-alpha*alpha*r2))/(sqrt(M_PI)*r) - erfc(alpha*r)/r2)/(2 * r)));
yher.push_back(hermite.Diff0(r2));
ynat.push_back(natural.Diff0(r2));
ynat2.push_back(natural2.Diff0(r2));
yquad.push_back(quadratic.Diff0(r2));
yhyp.push_back(hyperbolic.Diff0(r2));
ylin.push_back(linear.Diff0(r2));
ylin3.push_back(linear3.Diff0(r2));
ylin4.push_back(linear4.Diff0(r2));
xval.push_back(r2);
xsqrt.push_back(r);
}


Vector yex(yexact);
yex.PrintToFile("yexact.txt");
Vector yh(yher);
yh.PrintToFile("yhermite.txt");
Vector yn(ynat);
yn.PrintToFile("ynatural.txt");
Vector yn2(ynat2);
yn2.PrintToFile("ynatural2.txt");
Vector yq(yquad);
yq.PrintToFile("yquadratic.txt");
Vector yhy(yhyp);
yhy.PrintToFile("yhyperbolic.txt");
Vector yl(ylin);
yl.PrintToFile("ylinear.txt");
Vector yl3(ylin3);
yl3.PrintToFile("ylinear3.txt");
Vector yl4(ylin4);
yl4.PrintToFile("ylinear4.txt");
Vector xv(xval);
xv.PrintToFile("xval.txt");
Vector xs(xsqrt);
xs.PrintToFile("xsqrt.txt");

// mergetab xsqrt.txt:1  xval.txt:1  yexact.txt:1  yhermite.txt:1  yhyperbolic.txt:1  ylinear3.txt:1  ylinear4.txt:1  ylinear.txt:1  ynatural2.txt:1  ynatural.txt:1  yquadratic.txt:1 | plot -:2:"D-C" :"E-C" :"F-C" :"G-C" :"H-C" :"I-C" :"J-C"
return 0;
}
