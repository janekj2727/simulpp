/*
 *  Hermite cubic splines based on AbstractInterpolator interface.
 *  Originally for Simul++ simulation package
 *
 *  Implements the interface defined by AbstractInterpolator.
 *  Cubic splines with predefined values and first derivatives at all knots.
 *  For details see doxygen comments in HermiteCubSplines.hpp.
 *  Implementation must be here in header, because it is a template class.
 *  For inherited members must use this->, otherwise cannot compile...
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef HERMITECUBSPLINESHEADER
#define HERMITECUBSPLINESHEADER

#include "AbstractInterpolator.hpp"
#include <cmath>
#include <cassert>

/**
 * @class HermiteCubSplines
 * @ingroup interp
 * @brief Hermite cubic interpolation.
 *
 * This class provides Hermite cubic interpolation of given data.
 * Once initialized by function values and its first derivatives at grid points (knots),
 * the interpolator creates piecewise cubic polynomial to calculate function value in arbitrary point by @ref Eval() (or even @ref Eval0() when @ref minValue = 0.0).
 * The interpolation between two knots (\f$ x \in [x_i, x_{i+1}] \f$) is given by
 * \f{aligned}{
 *     t &= \left(x - h\lfloor x/h \rfloor\right)/h \\
 *     y(x) &= A_i + B_i t + C_i t^2 + D_i t^3
 * \f}
 * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
 *
 * The derivative estimation provided by @ref Diff() and @ref Diff0() is based on the formula
 * \f{equation}{
 *     y'(x)=\frac{1}{h}\left(B_i + 2 C_i t + 3 D_i t^2\right)
 * \f}
 * which is simply the derivative of the interpolating cubic polynomial.
 *
 * Parameters of cubic polynomials (\f$A_i\f$, \f$B_i\f$,...) are saved in the @ref params array
 * (which thus has the length of \f$4\times(\texttt{number_of_knots} - 1)\f$).
 * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ by
 * \f{aligned}{
 *     A_i &= f(x_i) \\
 *     B_i &= f'(x_i) \\
 *     C_i &= - 3 f(x_i) + 3 f(x_{i+1}) - 2 h f'(x_i) - h f'(x_{i+1}) \\
 *     D_i &= 2 f(x_i) - 2 f(x_{i+1}) + h f'(x_i) + h f'(x_{i+1})
 * \f}
 * This ensures that both function values and first derivatives at knots are equal to given values.
 * Resulting piecewise polynomial is continuous and has continuous first derivatives, second derivatives are not continuous.
 * Formulas inspired by MACSIMUS code @cite MACSIMUS_manual.
 * The interpolator can be defined with `double` or `long double` (or even `float`) precision.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */

/*
 * old params calculation
 * \f{aligned}{
 *     A_i =& \frac{1}{h^3} \left[(h - 2 x_i) x_{i+1}^2 f(x_{i}) + x_i^2 (3 h + 2 x_i) f(x_{i+1}) - x_i h x_{i+1}^2 f'(x_i) - x_{i+1} h x_i^2 f'(x_{i+1}) \right] \\
 *     B_i =& \frac{1}{h^3} \left[ 6 x_i x_{i+1} (f(_ix_{i}) - f(x_{i+1})) + h (3 x_i + h) x_{i+1} f'(x_i) + h x_i (3 x_i + 2 h) f'(x_{i+1})\right] \\
 *     C_i =& \frac{1}{h^3}\left[-3 (h + 2 x_i) (f(x_{i}) - f(x_{i+1})) - h (3 x_i + 2 h) f'(x_i) - h (3 x_i + h) f'(x_{i+1})\right] \\
 *     D_i =& \frac{1}{h^3} \left[2 (f(x_{i}) - f(x_{i+1})) + h (f'(x_i) + f'(x_{i+1}))\right]
 * \f}
 */
template <typename T>
class HermiteCubSplines : public AbstractInterpolator<T>
{
private:
    /**
     * @brief Default constructor of Hermite Cubic Splines object.
     *
     * Made private to disable its use.
     */
    HermiteCubSplines(){};

public:
    /**
     * @brief Destroy the Hermite Cubic Splines object
     *
     * Default and only destructor. Call to parent destructor should be enough.
     */
    ~HermiteCubSplines(){};
    /**
     * @brief Construct a new Hermite Cubic Splines object
     *
     * From given function values at @p N knots ( @p knotvalues) and values of first derivative at knots ( @p derivatives)
     * compute parameters of cubic polynomials and save them in @ref params.
     * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ by
     * \f{aligned}{
     *     A_i &= f(x_i) \\
     *     B_i &= f'(x_i) \\
     *     C_i &= - 3 f(x_i) + 3 f(x_{i+1}) - 2 h f'(x_i) - h f'(x_{i+1}) \\
     *     D_i &= 2 f(x_i) - 2 f(x_{i+1}) + h f'(x_i) + h f'(x_{i+1})
     * \f}
     * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param N Number of knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivatives at equidistant knots.
     */
    HermiteCubSplines(T *knotvalues, int N, T xMax, T xMin = 0.0L, T *derivatives = nullptr)
    {
        int i;
        this->size = N;
        this->params = new T[this->size * 4];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);

        for (i = 0; i < this->size - 1; i++)
        {
            this->params[4 * i + 0] = knotvalues[i];
            this->params[4 * i + 1] = this->step * derivatives[i];
            this->params[4 * i + 2] = -2.0L * this->step * derivatives[i] - 3.0L * knotvalues[i] + 3.0L * knotvalues[i + 1] - this->step * derivatives[i + 1];
            this->params[4 * i + 3] = 2.0L * knotvalues[i] + this->step * derivatives[i] - 2.0L * knotvalues[i + 1] + this->step * derivatives[i + 1];
        }
        // last point
        this->params[this->size * 4 - 4] = knotvalues[this->size - 1];
        this->params[this->size * 4 - 3] = 0.0;
        this->params[this->size * 4 - 2] = 0.0;
        this->params[this->size * 4 - 1] = 0.0;
    };
    /**
     * @brief Construct a new Hermite Cubic Splines object
     *
     * From given function values at knots ( @p knotvalues) and values of first derivative at knots ( @p derivatives)
     * compute parameters of cubic polynomials and save them in @ref params.
     * These parameters are computed for each interval between subsequent knot values \f$x_i\f$ and \f$x_{i+1}\f$ by
     * \f{aligned}{
     *     A_i &= f(x_i) \\
     *     B_i &= f'(x_i) \\
     *     C_i &= - 3 f(x_i) + 3 f(x_{i+1}) - 2 h f'(x_i) - h f'(x_{i+1}) \\
     *     D_i &= 2 f(x_i) - 2 f(x_{i+1}) + h f'(x_i) + h f'(x_{i+1})
     * \f}
     * where \f$h\f$ is the distance between the two knots \f$h = x_{i+1} - x_i\f$.
     * The interpolation range is given by @p xMin and @p xMax; both must be given, since @p derivatives must be defined for this type of interpolation.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives First derivatives at equidistant knots.
     */
    HermiteCubSplines(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
    {
        assert(knotvalues.size() == derivatives.size());
        int i;
        this->size = knotvalues.size();
        this->params = new T[this->size * 4];
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
        for (i = 0; i < this->size - 1; i++)
        {
            this->params[4 * i + 0] = knotvalues[i];
            this->params[4 * i + 1] = this->step * derivatives[i];
            this->params[4 * i + 2] = -2.0L * this->step * derivatives[i] - 3.0L * knotvalues[i] + 3.0L * knotvalues[i + 1] - this->step * derivatives[i + 1];
            this->params[4 * i + 3] = 2.0L * knotvalues[i] + this->step * derivatives[i] - 2.0L * knotvalues[i + 1] + this->step * derivatives[i + 1];
            // xi = this->minValue + i * this->step;
            // xip1 = xi + this->step;
            // this->params[4 * i + 0] = step_m3 * 2.0L * xi * (xi * xi * knotvalues[i + 1] - xip1 * xip1 * knotvalues[i]) + step_m2 * (xip1 * xip1 * knotvalues[i] + 3.0L * xi * xi * knotvalues[i + 1] - xi * xip1 * xip1 * derivatives[i] - xi * xi * xip1 * derivatives[i + 1]); // A_i
            // this->params[4 * i + 1] = step_m3 * (6.0L * xi * xip1 * (knotvalues[i] - knotvalues[i + 1])) + step_m2 * (3.0L * xi * (xip1 * derivatives[i] + xi * derivatives[i + 1])) + step_m1 * (xip1 * derivatives[i] + 2.0L * xi * derivatives[i + 1]);                        // B_i
            // this->params[4 * i + 2] = step_m3 * (6.0L * xi * (knotvalues[i + 1] - knotvalues[i])) - step_m2 * 3.0L * (knotvalues[i] - knotvalues[i + 1] + xi * (derivatives[i] + derivatives[i + 1])) - step_m1 * (2.0L * derivatives[i] + derivatives[i + 1]);                   // C_i
            // this->params[4 * i + 3] = step_m3 * (2.0L * (knotvalues[i] - knotvalues[i + 1])) + step_m2 * (derivatives[i] + derivatives[i + 1]); // D_i
        }
        // last point
        this->params[this->size * 4 - 4] = knotvalues[this->size - 1];
        this->params[this->size * 4 - 3] = 0.0;
        this->params[this->size * 4 - 2] = 0.0;
        this->params[this->size * 4 - 1] = 0.0;
    };
    /**
     * @brief Construct a new Hermite Cubic Splines object.
     *
     * Copy constructor for HermiteCubSplines.
     *
     * @param C Original interpolator that should be copied.
     */
    HermiteCubSplines(const HermiteCubSplines &C)
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
     * Assign the value of @p C to the left-hand-side instance of HermiteCubSplines.
     *
     * @note This is not a constructor!
     *
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    HermiteCubSplines &operator=(const HermiteCubSplines &C)
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
     * Evaluate the (approximation of) function value at the point @p xValue \f$= x \in [x_i, x_{i+1}]\f$.
     * Hermite piecewise cubic polynoms are used
     * \f{aligned}{
     *     t &= \left(x - h\lfloor x/h \rfloor\right)/h \\
     *     y(x) &= A_i + B_i t + C_i t^2 + D_i t^3
     * \f}
     * where \f$h\f$ is @ref step.
     * Parameters \f$A_i\f$, \f$B_i\f$,... are stored in @ref params.
     *
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    T Eval(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        xValue -= this->minValue;
        int i = (int)floor(xValue / this->step);
        xValue -= i * this->step;
        xValue /= this->step;
        i *= 4;
        return this->params[i] + xValue * (this->params[i + 1] + xValue * (this->params[i + 2] + xValue * this->params[i + 3]));
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue \f$= x \in [x_i, x_{i+1}]\f$.
     * This version assumes that the @ref minValue = 0.0.
     * Hermite piecewise cubic polynoms are used
     * \f{aligned}{
     *     t &= \left(x - h\lfloor x/h \rfloor\right)/h \\
     *     y(x) &= A_i + B_i t + C_i t^2 + D_i t^3
     * \f}
     * where \f$h\f$ is @ref step.
     * Parameters \f$A_i\f$, \f$B_i\f$,... are stored in @ref params.
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
        int i = (int)floor(xValue / this->step);
        xValue -= i * this->step;
        xValue /= this->step;
        i *= 4;
        return this->params[i] + xValue * (this->params[i + 1] + xValue * (this->params[i + 2] + xValue * this->params[i + 3]));
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * The following quadratic polynomial is used to approximate the derivative:
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=\frac{1}{h}\left(B_i + 2 C_i t + 3 D_i t^2\right)
     * \f}
     * where \f$h\f$ is the @ref step and \f$t=(x-x_i)/h\f$.
     * The formula can be derived by differentiating the interpolation cubic polynom.
     *
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = (int)floor((xValue - this->minValue) / this->step);
        xValue -= i * this->step + this->minValue;
        xValue /= this->step;
        i *= 4;
        return (this->params[i + 1] + xValue * (2.0L * this->params[i + 2] + 3.0L * xValue * this->params[i + 3])) / this->step;
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * The following quadratic polynomial is used to approximate the derivative:
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=\frac{1}{h}\left(B_i + 2 C_i t + 3 D_i t^2\right)
     * \f}
     * where \f$h\f$ is the @ref step and \f$t=(x-x_i)/h\f$.
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
        int i = (int)floor(xValue / this->step);
        xValue -= i * this->step;
        xValue /= this->step;
        i *= 4;
        return (this->params[i + 1] + xValue * (2.0L * this->params[i + 2] + 3.0L * xValue * this->params[i + 3])) / this->step;
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
        stream << indent << "- interpolation type: Hermite splines (see `HermiteCubSplines.hpp`)\n";
        stream << indent << "- usable interval: [" << this->minValue << ", " << this->maxValue << "]\n";
        stream << indent << "- resolution: " << (int) round(1.0/this->step) << " points per unit in argument (squared distance)\n";
    };
};

#endif
