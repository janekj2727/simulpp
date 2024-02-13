/*
 *  Natural cubic splines based on AbstractInterpolator interface.
 *  Originally for Simul++ simulation package
 *
 *  Implements the interface defined by AbstractInterpolator.
 *  Cubic splines with predefined values and continuous first and second derivatives.
 *  For details see doxygen comments in NaturalCubSplines.hpp.
 *  Implementation must be here in header, because it is a template class.
 *  For inherited members must use this->, otherwise cannot compile...
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef NATURALCUBSPLINESHEADER
#define NATURALCUBSPLINESHEADER

#include "AbstractInterpolator.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "LinearSystem.hpp"
#include <cmath>
#include <cassert>

/**
 * @class NaturalCubSplines
 * @ingroup interp
 * @brief Natural cubic interpolation.
 *
 * This class provides natural cubic splines interpolation of given data.
 * Once initialized by function values at grid points (knots) and first derivatives at initial and final knot,
 * the interpolator creates piecewise cubic polynomial to calculate function value in arbitrary point by @ref Eval() (or even @ref Eval0() when @ref minValue = 0.0).
 * The interpolation between two knots (\f$ x \in [x_i, x_{i+1}] \f$) is given by
 * \f{equation}{
 *     y(x) = A_i + B_i x + C_i x^2 + D_i x^3
 * \f}
 * where parameters of cubic polynomials (\f$A_i\f$, \f$B_i\f$,...) are saved in the @ref params array
 * (which thus has the length of \f$4\times(\texttt{number_of_knots} - 1)\f$).
 * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$
 * by solving global set of linear equations given by:
 * 1. Equality of interpolating polynomial and given function values at both knots
 * \f{aligned}{
 *     A_i + B_i x_i + C_i x_i^2 + D_i x_i^3 =& f(x_i) \\
 *     A_i + B_i x_{i+1} + C_i x_{i+1}^2 + D_i x_{i+1}^3 =& f(x_{i+1}) \quad \mathrm{for}~i = 0, 1,\dots, N-2
 * \f}
 * 2. Equality of first derivatives at knots (except for initial and final knots)
 * \f{equation}{
 *     B_i + 2 C_i x_{i+1} + 3 D_i x_{i+1}^2 - B_{i+1} - 2 C_{i+1} x_{i+1} - 3 D_{i+1} x_{i+1}^2 = 0 \quad \mathrm{for}~i = 0, 1,\dots, N-3
 * \f}
 * 3. Equality of second derivatives at knots (except for initial and final knots)
 * \f{equation}{
 *     2 C_i + 6 D_i x_{i+1} - 2 C_{i+1} - 6 D_{i+1} x_{i+1} = 0 \quad \mathrm{for}~i = 0, 1,\dots, N-3
 * \f}
 * 4. First derivatives at initial and final knots
 * \f{aligned}{
 *     B_0 + 2 C_0 x_0 + 3 D_0 x_0^2 =& f'(x_0) \\
 *     B_{N-1} + 2 C_{N-1} x_N + 3 D_{N-1} x_N^2 =& f'(x_N)
 * \f}
 * If the first derivatives are not given, the following condition applies:
 * Second derivatives at initial and final knots are equal to zero
 * \f{aligned}{
 *     2 C_0 + 6 D_0 x_0 =& 0 \\
 *     2 C_{N-1} + 6 D_{N-1} x_N =& 0
 * \f}
 *
 * The resulting linear system is solved by Gauss–Jordan elimination with pivoting (@ref LinearSystem::SolveGJ()).
 * Meanwhile, LinearSystem limited to `double` precision numbers.
 * Therefore, `long double` can be used but parameters are always computed with `double` precision.
 *
 * This construction ensures that function values are equal to defined values at knots and first and second derivatives are continuous at knots.
 * Moreover, first derivatives at the initial and final knots are equal to given values.
 *
 * The derivative estimation provided by @ref Diff() and @ref Diff0() is based on the formula
 * \f{equation}{
 *     f'(x)=B_i + 2 C_i x + 3 D_i x^2
 * \f}
 * which is simply the derivative of the interpolating cubic polynomial.
 *
 * The interpolator can be defined with `double` or `long double` (or even `float`) precision,
 * but precision of `long double` version is limited by the fact that LinearSystem can solve equations for parameters only in `double` precision.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */
template <typename T>
class NaturalCubSplines : public AbstractInterpolator<T>
{
private:
    /**
     * @brief Default constructor of Natural Cubic Splines object.
     *
     * Made private to disable its use.
     */
    NaturalCubSplines(){};

public:
    /**
     * @brief Destroy the Natural Cubic Splines object
     *
     * Default and only destructor. Call to parent destructor should be enough.
     */
    ~NaturalCubSplines(){};
    /**
     * @brief Construct a new Natural Cubic Splines object
     *
     * From given function values at @p N knots ( @p knotvalues) and values of first derivative at initial and final knot ( @p derivatives)
     * create system of linear equations for polynomial parameters,
     * solve the system and save parameters in @ref params.
     * If @p derivatives are not given, second derivatives at initial and final knots are set to 0.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param N Number of knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivatives at the initial and final knot.
     */
    NaturalCubSplines(T *knotvalues, int N, T xMax, T xMin = 0.0L, T *derivatives = nullptr)
    {
        int i, j;
        double xi;
        Matrix *A;
        Vector *b, *par;
        LinearSystem *linsystem;
        this->size = N;
        this->params = new T[this->size * 4];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        // construct system of linear equations A .  params = b
        A = new Matrix((this->size - 1) * 4, (this->size - 1) * 4);
        b = new Vector((this->size - 1) * 4);
        par = new Vector((this->size - 1) * 4);
        j = 0; // current row of A

        // the order of equations determines the size of round-off errors!!!
        /* ad 4. First derivatives at initial and final knots (if given)
         *     B_0 + 2 C_0 x_0 + 3 D_0 x_0^2 = f'(x_0)
         *     B_{N-1} + 2 C_{N-1} x_N + 3 D_{N-1} x_N^2 = f'(x_N)
         */
        xi = this->minValue; // x_i = x_0
        if (derivatives != nullptr)
        {
            A->operator()(j, 1) = 1.0;           // B_0
            A->operator()(j, 2) = 2.0 * xi;      // + 2 C_0 x_0
            A->operator()(j, 3) = 3.0 * xi * xi; // + 3 D_0 x_0^2
            b->operator()(j) = derivatives[0];   // = f'(x_0)
        }
        else
        /* ad 4 (if first derivatives not given):
         *       Second derivatives equal to zero.
         *       2 C_0 + 6 D_0 x_0 = 0
         *       2 C_{N-1} + 6 D_{N-1} x_N = 0
         */
        {
            A->operator()(j, 2) = 2.0;      // 2 C_0
            A->operator()(j, 3) = 6.0 * xi; // + 6 D_0 x_0
            b->operator()(j) = 0.0;         // = 0
        }
        j++; // next equation

        for (i = 0; i < this->size - 2; i++)
        {
            /* ad 1. Equality of interpolating polynomial and given function values at both knots
             *     A_i + B_i x_i + C_i x_i^2 + D_i x_i^3 = f(x_i)
             *     A_i + B_i x_{i+1} + C_i x_{i+1}^2 + D_i x_{i+1}^3 = f(x_{i+1})
             * for i = 0, 1,..., N-2
             */
            xi = this->minValue + i * this->step;
            A->operator()(j, i * 4 + 0) = 1.0;          // A_i
            A->operator()(j, i * 4 + 1) = xi;           // + B_i x_i
            A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_i^2
            A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_i^3
            b->operator()(j) = knotvalues[i];           // = f(x_i)
            xi += this->step;                           // x_i -> x_{i+1}
            j++;                                        // next equation
            A->operator()(j, i * 4 + 0) = 1.0;          // A_i
            A->operator()(j, i * 4 + 1) = xi;           // + B_i x_{i+1}
            A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_{i+1}^2
            A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_{i+1}^3
            b->operator()(j) = knotvalues[i + 1];       // = f(x_{i+1})
            j++;                                        // next equation

            /* ad 2. Equality of first derivatives at knots (except for initial and final knots)
             *     B_i + 2 C_i x_{i+1} + 3 D_i x_{i+1}^2 - B_{i+1} - 2 C_{i+1} x_{i+1} - 3 D_{i+1} x_{i+1}^2 = 0 i = 0, 1,..., N-3
             */
            A->operator()(j, i * 4 + 1) = 1.0;                  // B_i
            A->operator()(j, i * 4 + 2) = 2.0 * xi;             // + 2 C_i x_{i+1}
            A->operator()(j, i * 4 + 3) = 3.0 * xi * xi;        // + 3 D_i x_{i+1}^2
            A->operator()(j, (i + 1) * 4 + 1) = -1.0;           // -B_{i+1}
            A->operator()(j, (i + 1) * 4 + 2) = -2.0 * xi;      // - 2 C_{i+1} x_{i+1}
            A->operator()(j, (i + 1) * 4 + 3) = -3.0 * xi * xi; // - 3 D_{i+1} x_{i+1}^2
            b->operator()(j) = 0.0;                             // = 0
            j++;                                                // next equation

            /* ad3. Equality of second derivatives at knots (except for initial and final knots)
             *     2 C_i + 6 D_i x_{i+1} - 2 C_{i+1} - 6 D_{i+1} x_{i+1} = 0 for i = 0, 1,..., N-3
             */
            A->operator()(j, i * 4 + 2) = 2.0;             // 2 C_i
            A->operator()(j, i * 4 + 3) = 6.0 * xi;        // + 6 D_i x_{i+1}
            A->operator()(j, (i + 1) * 4 + 2) = -2.0;      // - 2 C_{i+1}
            A->operator()(j, (i + 1) * 4 + 3) = -6.0 * xi; // - 6 D_{i+1} x_{i+1}
            b->operator()(j) = 0.0;                        // = 0
            j++;                                           // next equation
        }

        i = this->size - 2;
        /* ad 1. Equality of interpolating polynomial and given function values at both knots
         *     A_i + B_i x_i + C_i x_i^2 + D_i x_i^3 = f(x_i)
         *     A_i + B_i x_{i+1} + C_i x_{i+1}^2 + D_i x_{i+1}^3 = f(x_{i+1})
         * for i = 0, 1,..., N-2
         */
        xi = this->minValue + i * this->step;
        A->operator()(j, i * 4 + 0) = 1.0;          // A_i
        A->operator()(j, i * 4 + 1) = xi;           // + B_i x_i
        A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_i^2
        A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_i^3
        b->operator()(j) = knotvalues[i];           // = f(x_i)
        xi += this->step;                           // x_i -> x_{i+1}
        j++;                                        // next equation
        A->operator()(j, i * 4 + 0) = 1.0;          // A_i
        A->operator()(j, i * 4 + 1) = xi;           // + B_i x_{i+1}
        A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_{i+1}^2
        A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_{i+1}^3
        b->operator()(j) = knotvalues[i + 1];       // = f(x_{i+1})
        j++;                                        // next equation

        /* ad 4. First derivatives at initial and final knots
         *     B_0 + 2 C_0 x_0 + 3 D_0 x_0^2 = f'(x_0)
         *     B_{N-1} + 2 C_{N-1} x_N + 3 D_{N-1} x_N^2 = f'(x_N)
         */
        xi = this->maxValue; // x_i = x_N
        if (derivatives != nullptr)
        {
            A->operator()(j, (this->size - 2) * 4 + 1) = 1.0;           // B_{N-1}
            A->operator()(j, (this->size - 2) * 4 + 2) = 2.0 * xi;      // + 2 C_{N-1} x_N
            A->operator()(j, (this->size - 2) * 4 + 3) = 3.0 * xi * xi; // + 3 D_{N-1} x_N^2
            b->operator()(j) = derivatives[1];                          // = f'(x_N)
        }
        else
        /* ad 4 (if first derivatives not given):
         *       Second derivatives equal to zero.
         *       2 C_0 + 6 D_0 x_0 = 0
         *       2 C_{N-1} + 6 D_{N-1} x_N = 0
         */
        {
            A->operator()(j, (this->size - 2) * 4 + 2) = 2.0;      // 2 C_{N-1}
            A->operator()(j, (this->size - 2) * 4 + 3) = 6.0 * xi; // + 6 D_{N-1} x_N
            b->operator()(j) = 0.0;                                // = 0
        }
        assert(j == A->GetNumberOfRows() - 1); // check if the number of equations equals the number of parameters

        // solve the system
        linsystem = new LinearSystem(*A, *b);
        *par = linsystem->SolveGJ();
        // assign values of parameters
        for (i = 0; i < (this->size - 1) * 4; i++)
        {
            this->params[i] = (T)par->Read(i);
        }
        // last point
        this->params[this->size * 4 - 4] = knotvalues[this->size - 1];
        this->params[this->size * 4 - 3] = 0.0;
        this->params[this->size * 4 - 2] = 0.0;
        this->params[this->size * 4 - 1] = 0.0;
        // delete auxiliary resources
        delete par;
        delete A;
        delete b;
        delete linsystem;
    };
    /**
     * @brief Construct a new Natural Cubic Splines object
     *
     * From given function values at knots ( @p knotvalues) and values of first derivative at knots ( @p derivatives)
     * create system of linear equations for polynomial parameters,
     * solve the system and save parameters in @ref params.
     * If @p derivatives are not given, second derivatives at initial and final knots are set to 0.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivatives at the initial and final knot.
     */
    NaturalCubSplines(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
    {
        int i, j;
        double xi;
        Matrix *A;
        Vector *b, *par;
        LinearSystem *linsystem;
        this->size = knotvalues.size();
        this->params = new T[this->size * 4];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        // construct system of linear equations A .  params = b
        A = new Matrix((this->size - 1) * 4, (this->size - 1) * 4);
        b = new Vector((this->size - 1) * 4);
        par = new Vector((this->size - 1) * 4);
        j = 0; // current row of A

        // the order of equations determines the size of round-off errors!!!
        /* ad 4. First derivatives at initial and final knots
         *     B_0 + 2 C_0 x_0 + 3 D_0 x_0^2 = f'(x_0)
         *     B_{N-1} + 2 C_{N-1} x_N + 3 D_{N-1} x_N^2 = f'(x_N)
         */
        xi = this->minValue; // x_i = x_0
        if (derivatives.size() > 0)
        {
            A->operator()(j, 1) = 1.0;           // B_0
            A->operator()(j, 2) = 2.0 * xi;      // + 2 C_0 x_0
            A->operator()(j, 3) = 3.0 * xi * xi; // + 3 D_0 x_0^2
            b->operator()(j) = derivatives[0];   // = f'(x_0)
        }
        else
        /* ad 4 (if first derivatives not given):
         *       Second derivatives equal to zero.
         *       2 C_0 + 6 D_0 x_0 = 0
         *       2 C_{N-1} + 6 D_{N-1} x_N = 0
         */
        {
            A->operator()(j, 2) = 2.0;      // 2 C_0
            A->operator()(j, 3) = 6.0 * xi; // + 6 D_0 x_0
            b->operator()(j) = 0.0;         // = 0
        }
        j++; // next equation

        for (i = 0; i < this->size - 2; i++)
        {
            /* ad 1. Equality of interpolating polynomial and given function values at both knots
             *     A_i + B_i x_i + C_i x_i^2 + D_i x_i^3 = f(x_i)
             *     A_i + B_i x_{i+1} + C_i x_{i+1}^2 + D_i x_{i+1}^3 = f(x_{i+1})
             * for i = 0, 1,..., N-2
             */
            xi = this->minValue + i * this->step;
            A->operator()(j, i * 4 + 0) = 1.0;          // A_i
            A->operator()(j, i * 4 + 1) = xi;           // + B_i x_i
            A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_i^2
            A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_i^3
            b->operator()(j) = knotvalues[i];           // = f(x_i)
            xi += this->step;                           // x_i -> x_{i+1}
            j++;                                        // next equation
            A->operator()(j, i * 4 + 0) = 1.0;          // A_i
            A->operator()(j, i * 4 + 1) = xi;           // + B_i x_{i+1}
            A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_{i+1}^2
            A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_{i+1}^3
            b->operator()(j) = knotvalues[i + 1];       // = f(x_{i+1})
            j++;                                        // next equation

            /* ad 2. Equality of first derivatives at knots (except for initial and final knots)
             *     B_i + 2 C_i x_{i+1} + 3 D_i x_{i+1}^2 - B_{i+1} - 2 C_{i+1} x_{i+1} - 3 D_{i+1} x_{i+1}^2 = 0 i = 0, 1,..., N-3
             */
            A->operator()(j, i * 4 + 1) = 1.0;                  // B_i
            A->operator()(j, i * 4 + 2) = 2.0 * xi;             // + 2 C_i x_{i+1}
            A->operator()(j, i * 4 + 3) = 3.0 * xi * xi;        // + 3 D_i x_{i+1}^2
            A->operator()(j, (i + 1) * 4 + 1) = -1.0;           // -B_{i+1}
            A->operator()(j, (i + 1) * 4 + 2) = -2.0 * xi;      // - 2 C_{i+1} x_{i+1}
            A->operator()(j, (i + 1) * 4 + 3) = -3.0 * xi * xi; // - 3 D_{i+1} x_{i+1}^2
            b->operator()(j) = 0.0;                             // = 0
            j++;                                                // next equation

            /* ad3. Equality of second derivatives at knots (except for initial and final knots)
             *     2 C_i + 6 D_i x_{i+1} - 2 C_{i+1} - 6 D_{i+1} x_{i+1} = 0 for i = 0, 1,..., N-3
             */
            A->operator()(j, i * 4 + 2) = 2.0;             // 2 C_i
            A->operator()(j, i * 4 + 3) = 6.0 * xi;        // + 6 D_i x_{i+1}
            A->operator()(j, (i + 1) * 4 + 2) = -2.0;      // - 2 C_{i+1}
            A->operator()(j, (i + 1) * 4 + 3) = -6.0 * xi; // - 6 D_{i+1} x_{i+1}
            b->operator()(j) = 0.0;                        // = 0
            j++;                                           // next equation
        }

        i = this->size - 2;
        /* ad 1. Equality of interpolating polynomial and given function values at both knots
         *     A_i + B_i x_i + C_i x_i^2 + D_i x_i^3 = f(x_i)
         *     A_i + B_i x_{i+1} + C_i x_{i+1}^2 + D_i x_{i+1}^3 = f(x_{i+1})
         * for i = 0, 1,..., N-2
         */
        xi = this->minValue + i * this->step;
        A->operator()(j, i * 4 + 0) = 1.0;          // A_i
        A->operator()(j, i * 4 + 1) = xi;           // + B_i x_i
        A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_i^2
        A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_i^3
        b->operator()(j) = knotvalues[i];           // = f(x_i)
        xi += this->step;                           // x_i -> x_{i+1}
        j++;                                        // next equation
        A->operator()(j, i * 4 + 0) = 1.0;          // A_i
        A->operator()(j, i * 4 + 1) = xi;           // + B_i x_{i+1}
        A->operator()(j, i * 4 + 2) = xi * xi;      // + C_i x_{i+1}^2
        A->operator()(j, i * 4 + 3) = xi * xi * xi; // + D_i x_{i+1}^3
        b->operator()(j) = knotvalues[i + 1];       // = f(x_{i+1})
        j++;                                        // next equation

        /* ad 4. First derivatives at initial and final knots
         *     B_0 + 2 C_0 x_0 + 3 D_0 x_0^2 = f'(x_0)
         *     B_{N-1} + 2 C_{N-1} x_N + 3 D_{N-1} x_N^2 = f'(x_N)
         */
        xi = this->maxValue; // x_i = x_N
        if (derivatives.size() > 0)
        {
            A->operator()(j, (this->size - 2) * 4 + 1) = 1.0;           // B_{N-1}
            A->operator()(j, (this->size - 2) * 4 + 2) = 2.0 * xi;      // + 2 C_{N-1} x_N
            A->operator()(j, (this->size - 2) * 4 + 3) = 3.0 * xi * xi; // + 3 D_{N-1} x_N^2
            b->operator()(j) = derivatives[1];                          // = f'(x_N)
        }
        else
        /* ad 4 (if first derivatives not given):
         *       Second derivatives equal to zero.
         *       2 C_0 + 6 D_0 x_0 = 0
         *       2 C_{N-1} + 6 D_{N-1} x_N = 0
         */
        {
            A->operator()(j, (this->size - 2) * 4 + 2) = 2.0;      // 2 C_{N-1}
            A->operator()(j, (this->size - 2) * 4 + 3) = 6.0 * xi; // + 6 D_{N-1} x_N
            b->operator()(j) = 0.0;                                // = 0
        }
        assert(j == A->GetNumberOfRows() - 1); // check if the number of equations equals the number of parameters

        // solve the system
        linsystem = new LinearSystem(*A, *b);
        *par = linsystem->SolveGJ();

        // assign values of parameters
        for (i = 0; i < (this->size - 1) * 4; i++)
        {
            this->params[i] = (T)par->Read(i);
        }
        // last point
        this->params[this->size * 4 - 4] = knotvalues[this->size - 1];
        this->params[this->size * 4 - 3] = 0.0;
        this->params[this->size * 4 - 2] = 0.0;
        this->params[this->size * 4 - 1] = 0.0;
        // delete auxiliary resources
        delete par;
        delete A;
        delete b;
        delete linsystem;
    };
    /**
     * @brief Construct a new Natural Cubic Splines object.
     *
     * Copy constructor for NaturalCubSplines.
     *
     * @param C Original interpolator that should be copied.
     */
    NaturalCubSplines(const NaturalCubSplines &C)
    {
        int i;
        this->size = C.size;
        this->params = new T[this->size * 4];
        for (i = 0; i < this->size * 4; i++)
        {
            this->params[i] = C.params[i];
        }
        this->minValue = C.minValue;
        this->maxValue = C.maxValue;
        this->step = C.step;
    };
    /**
     * @brief Default assignment operator.
     *
     * Assign the value of @p C to the left-hand-side instance of NaturalCubSplines.
     *
     * @note This is not a constructor!
     *
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    NaturalCubSplines &operator=(const NaturalCubSplines &C)
    {
        int i;
        assert(this->size >= C.size);
        this->size = C.size;
        for (i = 0; i < this->size * 4; i++)
        {
            this->params[i] = C.params[i];
        }
        this->minValue = C.minValue;
        this->maxValue = C.maxValue;
        this->step = C.step;
        return *this;
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue.
     * Hermite piecewise cubic polynoms are used
     * \f{equation}{
     *     y(x) = A_i + B_i x + C_i x^2 + D_i x^3 = A_i + x (B_i + x (C_i + x D_i))
     * \f}
     * where parameters \f$A_i\f$, \f$B_i\f$,... are stored in @ref params.
     *
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    T Eval(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = 4 * (int)floor((xValue - this->minValue) / this->step);
        return this->params[i] + xValue * (this->params[i + 1] + xValue * (this->params[i + 2] + xValue * this->params[i + 3]));
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * Hermite piecewise cubic polynoms are used
     * \f{equation}{
     *     y(x) = A_i + B_i x + C_i x^2 + D_i x^3 = A_i + x (B_i + x (C_i + x D_i))
     * \f}
     * where parameters \f$A_i\f$, \f$B_i\f$,... are stored in @ref params.
     *
     * @note No checks are provided that @ref minValue = 0.0 (only assertion in `DEBUG` version)!
     * Use to your own risk if you need to speed up the calculation a little.
     *
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    T Eval0(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        assert(this->minValue == 0.0);
        int i = 4 * (int)floor(xValue / this->step);
        return this->params[i] + xValue * (this->params[i + 1] + xValue * (this->params[i + 2] + xValue * this->params[i + 3]));
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * The following quadratic polynomial is used to approximate the derivative:
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=B_i + 2 C_i x + 3 D_i x^2
     * \f}
     * where \f$h\f$ is the @ref step.
     * The formula can be derived by differentiating the interpolation cubic polynom.
     *
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = 4 * (int)floor((xValue - this->minValue) / this->step);
        return this->params[i + 1] + xValue * (2.0L * this->params[i + 2] + 3.0L * xValue * this->params[i + 3]);
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * The following quadratic polynomial is used to approximate the derivative:
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=B_i + 2 C_i x + 3 D_i x^2
     * \f}
     * where \f$h\f$ is the @ref step.
     * The formula can be derived by differentiating the interpolation cubic polynom.
     *
     * @note No checks are provided that @ref minValue = 0.0 (only assertion in `DEBUG` version)!
     * Use to your own risk if you need to speed up the calculation a little.
     *
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff0(T xValue) const override
    {
        assert(this->minValue == 0.0);
        assert(xValue >= 0.0);
        assert(xValue <= this->maxValue);
        int i = 4 * (int)floor(xValue / this->step);
        return this->params[i + 1] + xValue * (2.0L * this->params[i + 2] + 3.0L * xValue * this->params[i + 3]);
    };
    /**
     * @brief Extend the range to zero ( @ref minValue = 0) to enable @ref Eval0().
     *
     * Extend the interpolation range to zero (set @ref minValue = 0.0) to enable Eval0().
     * The extension is possible only if @ref minValue is divisible by @ref step without remainder.
     * The first cubic polynomial of the original interpolation is used for extrapolation.
     *
     * @return Number of added knots (or equivalent number) or -1 if the extension cannot be done.
     */
    int ExtendToZero() override
    {
        int i;
        T *old_params = this->params;
        int no_added = (int)roundl(this->minValue / this->step);
        if ((T)no_added != this->minValue / this->step)
        {
            return -1;
        }
        this->size += no_added;
        this->params = new T[this->size * 4];
        for (i = this->size * 4 - 1; i >= no_added * 4; i--)
        {
            this->params[i] = old_params[i - no_added * 4];
        }
        for (i = no_added * 4 - 1; i >= 0; i--)
        {
            this->params[i] = this->params[i + 4];
        }
        this->minValue = 0.0;
        delete[] old_params;

        return no_added;
    };
    void PrintInfo(std::ofstream &stream, int indentation = 0) const override
    {
        std::string indent((size_t) indentation, ' ');
        stream << indent << "- interpolation type: natural cubic splines (see `NaturalCubSplines.hpp`)\n";
        stream << indent << "- usable interval: [" << this->minValue << ", " << this->maxValue << "]\n";
        stream << indent << "- resolution: " << (int) round(1.0/this->step) << " points per unit in argument (squared distance)\n";
    };
};

#endif
