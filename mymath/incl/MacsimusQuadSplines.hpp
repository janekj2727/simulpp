/*
 *  MACSIMUS-style quadratic splines based on AbstractInterpolator interface.
 *  Originally for Simul++ simulation package
 *
 *  Implements the interface defined by AbstractInterpolator.
 *  Quadratic splines with predefined values at all knots and first derivative at a final knot.
 *  For details see doxygen comments in MacsimusQuadSplines.hpp.
 *  Implementation must be here in header, because it is a template class.
 *  For inherited members must use this->, otherwise cannot compile...
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef MACSIMUSQUADSPLINESHEADER
#define MACSIMUSQUADSPLINESHEADER

#include "AbstractInterpolator.hpp"
#include <cmath>
#include <cassert>

/**
 * @class MacsimusQuadSplines
 * @ingroup interp
 * @brief MACSIMUS-style quadratic interpolation.
 *
 * This class provides a MACSIMUS-style quadratic interpolation of given data.
 * Once initialized by function values at grid points (knots) and the first derivative in the final knot,
 * the interpolator creates piecewise quadratic polynomial to calculate function value in arbitrary point by @ref Eval() (or even @ref Eval0() when @ref minValue = 0.0).
 * The interpolation between two knots (\f$ x \in [x_i, x_{i+1}] \f$) is given by
 * \f{equation}{
 *     y(x) = A_i + B_i x + C_i x^2
 * \f}
 * where parameters of quadratic polynomials (\f$A_i\f$, \f$B_i\f$,...) are saved in the @ref params array
 * (which thus has the length of \f$3\times(\texttt{number_of_knots} - 1)\f$).
 * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ by
 * \f{aligned}{
 *     C_i =& \left(hf'(x_{i+1}) + f(x_i) - f(x_{i+1}\right)/h^2 \\
 *     B_i =& f'(x_{i+1}) - 2 C_i x_{i+1} \\
 *     A_i =& f(x_i) - C_i x_i^2 - B_i x_i \\
 *     f'(x_i) =& 2C_i x_i + B_i
 * \f}
 * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
 * Parameter values are computed from the end of the interpolation interval.
 * Starting value of the first derivative at the last point is needed.
 * This computation ensures that both function values and first derivatives at knots are continuous, second derivatives are not continuous.
 * Function values at knots are equal to the given values.
 *
 * The derivative estimation provided by @ref Diff() and @ref Diff0() is based on the formula
 * \f{equation}{
 *     f'(x)=B_i + 2 C_i x
 * \f}
 * which is simply the derivative of the interpolating quadratic polynomial.
 *
 * Formulas inspired by MACSIMUS code @cite MACSIMUS_manual.
 * The interpolator can be defined with `double` or `long double` (or even `float`) precision.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */
template <typename T>
class MacsimusQuadSplines : public AbstractInterpolator<T>
{
private:
    /**
     * @brief Default constructor of MACSIMUS Quadratic Splines object.
     *
     * Made private to disable its use.
     */
    MacsimusQuadSplines(){};

public:
    /**
     * @brief Destroy the MACSIMUS Quadratic Splines object
     *
     * Default and only destructor. Call to parent destructor should be enough.
     */
    ~MacsimusQuadSplines(){};
    /**
     * @brief Construct a new MACSIMUS Quadratic Splines object
     *
     * From given function values at @p N knots ( @p knotvalues) and the value of first derivative at final knot ( @p derivatives),
     * compute parameters of quadratic polynomials and save them in @ref params.
     * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ starting from the last by
     * \f{aligned}{
     *     C_i =& \left(hf'(x_{i+1}) + f(x_i) - f(x_{i+1}\right)/h^2 \\
     *     B_i =& f'(x_{i+1}) - 2 C_i x_{i+1} \\
     *     A_i =& f(x_i) - C_i x_i^2 - B_i x_i \\
     *     f'(x_i) =& 2C_i x_i + B_i
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
    MacsimusQuadSplines(T *knotvalues, int N, T xMax, T xMin = 0.0L, T *derivatives = nullptr)
    {
        int i;
        T step_m2, xi, xip1, dfxi; // 1/h^2, x_i, x_{i+1}, f'(x_i)
        this->size = N;
        this->params = new T[this->size * 3];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        step_m2 = pow(this->step, -2.0);
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
            this->params[3 * i + 2] = step_m2 * (this->step * dfxi + knotvalues[i] - knotvalues[i + 1]);                // C_i
            this->params[3 * i + 1] = dfxi - 2.0 * this->params[3 * i + 2] * xip1;                                      // B_i
            this->params[3 * i + 0] = knotvalues[i] - this->params[3 * i + 2] * xi * xi - this->params[3 * i + 1] * xi; // A_i
            dfxi = 2.0 * this->params[3 * i + 2] * xi + this->params[3 * i + 1];
        }
    };
    /**
     * @brief Construct a new MACSIMUS Quadratic Splines object
     *
     * From given function values at knots ( @p knotvalues) and the value of first derivative at final knot ( @p derivatives),
     * compute parameters of quadratic polynomials and save them in @ref params.
     * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ starting from the last by
     * \f{aligned}{
     *     C_i =& \left(hf'(x_{i+1}) + f(x_i) - f(x_{i+1}\right)/h^2 \\
     *     B_i =& f'(x_{i+1}) - 2 C_i x_{i+1} \\
     *     A_i =& f(x_i) - C_i x_i^2 - B_i x_i \\
     *     f'(x_i) =& 2C_i x_i + B_i
     * \f}
     * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined (and contain exactly one value) for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivative at final knot.
     */
    MacsimusQuadSplines(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
    {
        int i;
        T step_m2, xi, xip1, dfxi; // 1/h^2, x_i, x_{i+1}, f'(x_i)
        this->size = knotvalues.size();
        this->params = new T[this->size * 3];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        step_m2 = pow(this->step, -2.0);
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
            this->params[3 * i + 2] = step_m2 * (this->step * dfxi + knotvalues[i] - knotvalues[i + 1]);                // C_i
            this->params[3 * i + 1] = dfxi - 2.0 * this->params[3 * i + 2] * xip1;                                      // B_i
            this->params[3 * i + 0] = knotvalues[i] - this->params[3 * i + 2] * xi * xi - this->params[3 * i + 1] * xi; // A_i
            dfxi = 2.0 * this->params[3 * i + 2] * xi + this->params[3 * i + 1];
        }
    };
    /**
     * @brief Construct a new MACSIMUS Quadratic Splines object.
     *
     * Copy constructor for MacsimusQuadSplines.
     *
     * @param C Original interpolator that should be copied.
     */
    MacsimusQuadSplines(const MacsimusQuadSplines &C)
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
     * Assign the value of @p C to the left-hand-side instance of MacsimusQuadSplines.
     *
     * @note This is not a constructor!
     *
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    MacsimusQuadSplines &operator=(const MacsimusQuadSplines &C)
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
     * Hermite piecewise quadratic polynoms are used
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
        int i = 3 * (int)floor((xValue - this->minValue) / this->step);
        return this->params[i] + xValue * (this->params[i + 1] + xValue * this->params[i + 2]);
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * Hermite piecewise quadratic polynoms are used
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
        int i = 3 * (int)floor(xValue / this->step);
        return this->params[i] + xValue * (this->params[i + 1] + xValue * this->params[i + 2]);
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * The following linear polynomial is used to approximate the derivative
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=B_i + 2 C_i x
     * \f}
     * where \f$h\f$ is the @ref step.
     * The formula can be derived by differentiating the interpolation quadratic polynom.
     *
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = 3 * (int)floor((xValue - this->minValue) / this->step);
        return this->params[i + 1] + xValue * 2.0L * this->params[i + 2];
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * The following linear polynomial is used to approximate the derivative
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=B_i + 2 C_i x
     * \f}
     * where \f$h\f$ is the @ref step.
     * The formula can be derived by differentiating the interpolation quadratic polynom.
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
        return this->params[i + 1] + xValue * 2.0L * this->params[i + 2];
    };
    /**
     * @brief Extend the range to zero ( @ref minValue = 0) to enable @ref Eval0().
     *
     * Extend the interpolation range to zero (set @ref minValue = 0.0) to enable Eval0().
     * The extension is possible only if @ref minValue is divisible by @ref step without remainder.
     * The first quadratic function of the original interpolation is used for extrapolation.
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
        stream << indent << "- interpolation type: quadratic splines (MACSIMUS type, see `MacsimusQuadSplines.hpp`)\n";
        stream << indent << "- usable interval: [" << this->minValue << ", " << this->maxValue << "]\n";
        stream << indent << "- resolution: " << (int) round(1.0/this->step) << " points per unit in argument (squared distance)\n";
    };
};

#endif
