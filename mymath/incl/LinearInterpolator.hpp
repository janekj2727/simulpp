/*
 *  Basic linear interpolation based on AbstractInterpolator interface.
 *  Originally for Simul++ simulation package
 *
 *  Implements the interface defined by AbstractInterpolator.
 *  For details see doxygen comments in LinearInterpolator.hpp.
 *  Implementation must be here in header, because it is a template class.
 *  For inherited members must use this->, otherwise cannot compile...
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef LINEARINTERPOLATORHEADER
#define LINEARINTERPOLATORHEADER

#include "AbstractInterpolator.hpp"
#include <cmath>
#include <cassert>

/**
 * @class LinearInterpolator
 * @ingroup interp
 * @brief Basic 2-point linear interpolation.
 *
 * This class provides the basic linear interpolation of given data.
 * Once initialized by function values at grid points (knots), the interpolator returns function value in arbitrary point by @ref Eval() (or even @ref Eval0() when @ref minValue = 0.0).
 * The interpolation between two knots (\f$ x \in [x_i, x_{i+1}] \f$) is given by
 * \f{equation}{
 *     y(x) = y_i + \frac{y_{i+1} - y_{i}}{x_{i+1} - x_i} (x - x_i) = y_i + \frac{y_{i+1} - y_{i}}{\mathtt{step}} (x - x_i)
 * \f}
 * The function values themselves are saved in the @ref params array.
 * 
 *  * The derivative estimation provided by @ref Diff() and @ref Diff0() is based on the simple formula
 * \f{equation}{
 *     f'(x)=\frac{y_{i+1} - y_{i}}{\mathtt{step}}
 * \f}
 * 
 * The interpolator can be defined with `double` or `long double` (or even `float`) precision.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */
template <typename T>
class LinearInterpolator : public AbstractInterpolator<T>
{
private:
    /**
     * @brief Default constructor of Linear Interpolator object.
     *
     * Made private to disable its use.
     */
    LinearInterpolator(){};

public:
    /**
     * @brief Destroy the Linear Interpolator object
     *
     * Default and only destructor. Call to parent destructor should be enough.
     */
    ~LinearInterpolator(){};
    /**
     * @brief Construct a new Linear Interpolator object
     *
     * Define function values at knots in @p knotvalues and the number of knots in @p N.
     * The interpolation range is given by @p xMin and @p xMax; the former is assumed to be 0.0 if not specified otherwise.
     * Derivatives provided by @p derivatives are ignored, even if defined.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param N Number of knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives Ignored.
     */
    LinearInterpolator(T *knotvalues, int N, T xMax, T xMin = 0.0L, T *derivatives = nullptr)
    {
        int i;
        this->size = N;
        this->params = new T[this->size];
        for (i = 0; i < this->size; i++)
        {
            this->params[i] = knotvalues[i];
        }
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
    };
    /**
     * @brief Construct a new Linear Interpolator object
     *
     * Define function values at knots in @p knotvalues; the number of knots can be determined by @p knotvalues size.
     * The interpolation range is given by @p xMin and @p xMax; the former is assumed to be 0.0 if not specified otherwise.
     * Function derivatives provided by @p derivatives are ignored.
     *
     * @param knotvalues Function values at equidistant knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives Ignored.
     */
    LinearInterpolator(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
    {
        int i;
        this->size = knotvalues.size();
        this->params = new T[this->size];
        for (i = 0; i < this->size; i++)
        {
            this->params[i] = knotvalues[i];
        }
        this->minValue = xMin;
        this->maxValue = xMax;
        this->step = (this->maxValue - this->minValue) / (this->size - 1);
    };
    /**
     * @brief Construct a new Linear Interpolator object.
     *
     * Copy constructor for LinearInterpolator.
     *
     * @param C Original interpolator that should be copied.
     */
    LinearInterpolator(const LinearInterpolator &C)
    {
        int i;
        this->size = C.size;
        this->params = new T[this->size];
        for (i = 0; i < this->size; i++)
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
     * Assign the value of @p C to the left-hand-side instance of LinearInterpolator.
     *
     * @note This is not a constructor!
     *
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    LinearInterpolator &operator=(const LinearInterpolator &C)
    {
        int i;
        this->size = C.size;
        for (i = 0; i < this->size; i++)
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
     * Basic linear interpolation is used:
     * \f{equation}{
     *     y(x) = y_i + \frac{y_{i+1} - y_{i}}{x_{i+1} - x_i}(x - x_i) = y_i + \frac{y_{i+1} - y_{i}}{\mathtt{step}}(x - x_i)
     * \f}
     *
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    T Eval(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = (int)floor((xValue - this->minValue) / this->step);
        return this->params[i] + (this->params[i + 1] - this->params[i]) / this->step * (xValue - (this->minValue + i * this->step));
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * Basic linear interpolation is used:
     * \f{equation}{
     *     y(x) = y_i + \frac{y_{i+1} - y_{i}}{x_{i+1} - x_i} (x - x_i) = y_i + \frac{y_{i+1} - y_{i}}{\mathtt{step}}(x - i * \mathtt{step})
     * \f}
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
        return this->params[i] + (this->params[i + 1] - this->params[i]) / this->step * (xValue - i * this->step);
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     * 
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * The most basic formula to approximate the derivative is used: 
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=\frac{f(x_{i+1}) - f(x_i)}{h}
     * \f} 
     * where \f$h\f$ is the @ref step.
     * The resulting derivatives are not continuous (piecewise constant function).
     * 
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue);
        int i = (int)floor((xValue - this->minValue) / this->step);
        return (this->params[i + 1] - this->params[i]) / this->step;
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     * 
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0 (\f$x_0 = 0.0\f$).
     * The most basic formula to approximate the derivative is used: 
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=\frac{f(x_{i+1}) - f(x_i)}{h}
     * \f} 
     * where \f$h\f$ is the @ref step.
     * The resulting derivatives are not continuous (piecewise constant function).
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
        return (this->params[i + 1] - this->params[i]) / this->step;
    };
     /**
     * @brief Extend the range to zero ( @ref minValue = 0) to enable @ref Eval0().
     * 
     * Extend the interpolation range to zero (set @ref minValue = 0.0) to enable Eval0().
     * The extension is possible only if @ref minValue is divisible by @ref step without remainder.
     * Basic linear extrapolation from the first two points is used.
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
        this->params = new T[this->size];
        for (i = this->size - 1; i >= no_added; i--)
        {
            this->params[i] = old_params[i - no_added];
        }
        // extrapolation
        for (i = no_added - 1; i >= 0; i--)
        {
            this->params[i] = old_params[0] - (old_params[1] - old_params[0]) * (no_added - i);
        }
        this->minValue = 0.0;
        delete[] old_params;

        return no_added;
    };
    void PrintInfo(std::ofstream &stream, int indentation = 0) const override
    {
        std::string indent((size_t) indentation, ' ');
        stream << indent << "- interpolation type: basic linear interpolation (see `LinearInterpolator.hpp`)\n";
        stream << indent << "- usable interval: [" << this->minValue << ", " << this->maxValue << "]\n";
        stream << indent << "- resolution: " << (int) round(1.0/this->step) << " points per unit in argument (squared distance)\n";
    };
};

#endif
