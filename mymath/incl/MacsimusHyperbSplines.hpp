/*
 *  MACSIMUS-style hyperbolic splines based on AbstractInterpolator interface.
 *  Originally for Simul++ simulation package
 *
 *  Implements the interface defined by AbstractInterpolator.
 *  Hyperbolic splines with predefined values at all knots and first derivative at a final knot.
 *  For details see doxygen comments in MacsimusHyperbSplines.hpp.
 *  Implementation must be here in header, because it is a template class.
 *  For inherited members must use this->, otherwise cannot compile...
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef MACSIMUSHYPERBSPLINESHEADER
#define MACSIMUSHYPERBSPLINESHEADER

#include "AbstractInterpolator.hpp"
#include <cmath>
#include <cassert>

/**
 * @class MacsimusHyperbSplines
 * @ingroup interp
 * @brief MACSIMUS-style hyperbolic interpolation.
 *
 * This class provides a MACSIMUS-style hyperbolic interpolation of given data.
 * Once initialized by function values at grid points (knots) and the first derivative in the final knot,
 * the interpolator creates piecewise hyperbolic polynomial to calculate function value in arbitrary point by @ref Eval() (or even @ref Eval0() when @ref minValue = 0.0).
 * The interpolation between two knots (\f$ x \in [x_i, x_{i+1}] \f$) is given by
 * \f{equation}{
 *     y(x) = A_i + \frac{B_i}{C_i + x}
 * \f}
 * where parameters of hyperbolic polynomials (\f$A_i\f$, \f$B_i\f$,...) are saved in the @ref params array
 * (which thus has the length of \f$3\times(\texttt{number_of_knots}-1\f$).
 * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ by
 * \f{aligned}{
 *     C_i =& \frac{x_i(f(x_{i+1}) - f(x_i))- x_{i+1} h f'(x_{i+1})}{h f'(x_{i+1}) - (f(x_{i+1}) - f(x_i))} \\
 *     B_i =& -f'(x_{i+1}) (C + x_{i+1})^2 \\
 *     A_i =&  f(x_{i+1}) + f'(x_{i+1}) (C + x_{i+1})\\
 *     f'(x_i) =& -\frac{B}{(C+x_i)^2}
 * \f}
 * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
 * Parameter values are computed from the end of the interpolation interval.
 * Starting value of the first derivative at the last point is needed.
 * This computation ensures that both function values and first derivatives at knots are continuous, second derivatives are not continuous.
 * Function values at knots are equal to the given values.
 *
 * The derivative estimation provided by @ref Diff() and @ref Diff0() is based on the formula
 * \f{equation}{
 *     y'(x) = -\frac{B_i}{(C_i + x)^2}
 * \f}
 * which is simply the derivative of the interpolating hyperbolic function.
 *
 * @note This interpolation is suitable only for monotonous functions!
 *
 * Formulas inspired by MACSIMUS code @cite MACSIMUS_manual.
 * The interpolator can be defined with `double` or `long double` (or even `float`) precision.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */
template <typename T>
class MacsimusHyperbSplines : public AbstractInterpolator<T>
{
private:
    /**
     * @brief Default constructor of MACSIMUS Hyperbolic Splines object.
     *
     * Made private to disable its use.
     */
    MacsimusHyperbSplines(){};

public:
    /**
     * @brief Destroy the MACSIMUS Hyperbolic Splines object
     *
     * Default and only destructor. Call to parent destructor should be enough.
     */
    ~MacsimusHyperbSplines(){};
    /**
     * @brief Construct a new MACSIMUS Hyperbolic Splines object
     *
     * From given function values at @p N knots ( @p knotvalues) and the value of first derivative at final knot ( @p derivatives),
     * compute parameters of hyperbolic polynomials and save them in @ref params.
     * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ starting from the last by
     * \f{aligned}{
     *     C_i =& \frac{x_i(f(x_{i+1}) - f(x_i))- x_{i+1} h f'(x_{i+1})}{h f'(x_{i+1}) - (f(x_{i+1}) - f(x_i))} \\
     *     B_i =& -f'(x_{i+1}) (C + x_{i+1})^2 \\
     *     A_i =&  f(x_{i+1}) + f'(x_{i+1}) (C + x_{i+1})\\
     *     f'(x_i) =& -\frac{B}{(C+x_i)^2}
     * \f}
     * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined (and contain exactly one value) for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param N Number of knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivative at final knot.
     */
    MacsimusHyperbSplines(T *knotvalues, int N, T xMax, T xMin = 0.0L, T *derivatives = nullptr)
    {
        int i;
        T xi, xip1, dfxi; // x_i, x_{i+1}, f'(x_i)
        this->size = N;
        this->params = new T[this->size * 3];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        if (derivatives == nullptr)
        {
            dfxi = 0.0; // just to provide suitable value
        }
        else
        {
            dfxi = derivatives[0];
        }
        // last point
        this->params[this->size * 3 - 3] = knotvalues[this->size - 1];
        this->params[this->size * 3 - 2] = 0.0;
        this->params[this->size * 3 - 1] = 0.0;
        // from the last point to xMin
        for (i = this->size - 2; i >= 0; i--)
        {
            xi = this->minValue + i * this->step;
            xip1 = xi + this->step;
            this->params[3 * i + 2] = (xi * (knotvalues[i + 1] - knotvalues[i]) - xip1 * this->step * dfxi)(this->step * dfxi - (knotvalues[i + 1] - knotvalues[i])); // C_i
            this->params[3 * i + 1] = -dfxi * (this->params[3 * i + 2] + xip1) * (this->params[3 * i + 2] + xip1);                                                    // B_i
            this->params[3 * i + 0] = knotvalues[i + 1] + dfxi * (this->params[3 * i + 2] + xip1);                                                                    // A_i
            dfxi = -this->params[3 * i + 1] / (this->params[3 * i + 2] + xi) / (this->params[3 * i + 2] + xi);
        }
    };
    /**
     * @brief Construct a new MACSIMUS Hyperbolic Splines object
     *
     * From given function values at knots ( @p knotvalues) and the value of first derivative at final knot ( @p derivatives),
     * compute parameters of hyperbolic polynomials and save them in @ref params.
     * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ starting from the last by
     * \f{aligned}{
     *     C_i =& \frac{x_i(f(x_{i+1}) - f(x_i))- x_{i+1} h f'(x_{i+1})}{h f'(x_{i+1}) - (f(x_{i+1}) - f(x_i))} \\
     *     B_i =& -f'(x_{i+1}) (C + x_{i+1})^2 \\
     *     A_i =&  f(x_{i+1}) + f'(x_{i+1}) (C + x_{i+1})\\
     *     f'(x_i) =& -\frac{B}{(C+x_i)^2}
     * \f}
     * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined (and contain exactly one value) for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivative at final knot.
     */
    MacsimusHyperbSplines(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
    {
        int i;
        T xi, xip1, dfxi; // x_i, x_{i+1}, f'(x_i)
        this->size = knotvalues.size();
        this->params = new T[this->size * 3];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        if (derivatives.size() < 1)
        {
            dfxi = 0.0; // just to provide suitable value
        }
        else
        {
            dfxi = derivatives[0];
        }
        // last point
        this->params[this->size * 3 - 3] = knotvalues[this->size - 1];
        this->params[this->size * 3 - 2] = 0.0;
        this->params[this->size * 3 - 1] = 0.0;
        // from the last point to xMin
        for (i = this->size - 2; i >= 0; i--)
        {
            xi = this->minValue + i * this->step;
            xip1 = xi + this->step;
            this->params[3 * i + 2] = (xi * (knotvalues[i + 1] - knotvalues[i]) - xip1 * this->step * dfxi) / (this->step * dfxi - (knotvalues[i + 1] - knotvalues[i])); // C_i
            this->params[3 * i + 1] = -dfxi * (this->params[3 * i + 2] + xip1) * (this->params[3 * i + 2] + xip1);                                                       // B_i
            this->params[3 * i + 0] = knotvalues[i + 1] + dfxi * (this->params[3 * i + 2] + xip1);                                                                       // A_i
            dfxi = -this->params[3 * i + 1] / (this->params[3 * i + 2] + xi) / (this->params[3 * i + 2] + xi);
        }
    };
    /**
     * @brief Construct a new MACSIMUS Hyperbolic Splines object.
     *
     * Copy constructor for MacsimusHyperbSplines.
     *
     * @param C Original interpolator that should be copied.
     */
    MacsimusHyperbSplines(const MacsimusHyperbSplines &C)
    {
        int i;
        this->size = C.size;
        this->params = new T[this->size * 3];
        for (i = 0; i < this->size * 3; i++)
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
     * Assign the value of @p C to the left-hand-side instance of MacsimusHyperbSplines.
     *
     * @note This is not a constructor!
     *
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    MacsimusHyperbSplines &operator=(const MacsimusHyperbSplines &C)
    {
        int i;
        assert(this->size >= C.size);
        this->size = C.size;
        for (i = 0; i < this->size * 3; i++)
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
     * Hermite piecewise hyperbolic polynoms are used
     * \f{equation}{
     *     y(x) = A_i + \frac{B_i}{C_i + x}
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
        int i = 3 * (int)floor((xValue - this->minValue) / this->step);
        return this->params[i] + this->params[i + 1] / (xValue + this->params[i + 2]);
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * Hermite piecewise hyperbolic polynoms are used
     * \f{equation}{
     *     y(x) = A_i + \frac{B_i}{C_i + x}
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
        int i = 3 * (int)floor(xValue / this->step);
        return this->params[i] + this->params[i + 1] / (xValue + this->params[i + 2]);
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * The following formula is used to approximate the derivative
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     y'(x) = -\frac{B_i}{(C_i + x)^2}
     * \f}
     * where \f$h\f$ is the @ref step.
     * The formula can be derived by differentiating the interpolation hyperbolic function.
     *
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = 3 * (int)floor((xValue - this->minValue) / this->step);
        return -this->params[i + 1] / ((xValue + this->params[i + 2]) * (xValue + this->params[i + 2]));
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * The following formula is used to approximate the derivative
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     y'(x) = -\frac{B_i}{(C_i + x)^2}
     * \f}
     * where \f$h\f$ is the @ref step.
     * The formula can be derived by differentiating the interpolation hyperbolic function.
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
        int i = 3 * (int)floor(xValue / this->step);
        return -this->params[i + 1] / ((xValue + this->params[i + 2]) * (xValue + this->params[i + 2]));
    };
    /**
     * @brief Extend the range to zero ( @ref minValue = 0) to enable @ref Eval0().
     *
     * Extend the interpolation range to zero (set @ref minValue = 0.0) to enable Eval0().
     * The extension is possible only if @ref minValue is divisible by @ref step without remainder.
     * The first hyperbolic function of the original interpolation is used for extrapolation.
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
        this->params = new T[this->size * 3];
        for (i = this->size * 3 - 1; i >= no_added * 3; i--)
        {
            this->params[i] = old_params[i - no_added * 3];
        }
        for (i = no_added * 3 - 1; i >= 0; i--)
        {
            this->params[i] = this->params[i + 3];
        }
        this->minValue = 0.0L;
        delete[] old_params;

        return no_added;
    };
    void PrintInfo(std::ofstream &stream, int indentation = 0) const override
    {
        std::string indent((size_t) indentation, ' ');
        stream << indent << "- interpolation type: hyperbolic splines (MACSIMUS type, see `MacsimusHyperbSplines.hpp`)\n";
        stream << indent << "- usable interval: [" << this->minValue << ", " << this->maxValue << "]\n";
        stream << indent << "- resolution: " << (int) round(1.0/this->step) << " points per unit in argument (squared distance)\n";
    };
};

#endif
