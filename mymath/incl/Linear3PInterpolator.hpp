/*
 *  3-point linear interpolation (DL_POLY style) based on AbstractInterpolator interface.
 *  Originally for Simul++ simulation package
 *
 *  Implements the interface defined by AbstractInterpolator.
 *  For details see doxygen comments in Linear3PInterpolator.hpp.
 *  Implementation must be here in header, because it is a template class.
 *  For inherited members must use this->, otherwise cannot compile...
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef LINEAR3PINTERPOLATORHEADER
#define LINEAR3PINTERPOLATORHEADER

#include "AbstractInterpolator.hpp"
#include <cmath>
#include <cassert>

/**
 * @class Linear3PInterpolator
 * @ingroup interp
 * @brief Three-point linear interpolation.
 *
 * This class provides the 3-point linear interpolation of given data.
 * This type of linear interpolation is based on DL_POLY manual @cite dl_poly4-manual and code.
 * Once initialized by function values at grid points (knots), the interpolator returns function value in arbitrary point by @ref Eval() (or even @ref Eval0() when @ref minValue = 0.0).
 * The interpolation between two knots (\f$ x \in [x_i, x_{i+1}] \f$) uses three consequent knots (\f$y_i\f$, \f$y_{i+1}\f$ and \f$y_{i+2}\f$) and is given by
 * \f{aligned}{
 *   i =& \lfloor (x - x_0) / h \rfloor \\
 *   t =& \frac{x - x_i}{h} \\
 *   s_1 =& y_i + (y_{i+1} - y_i)t \\
 *   s_2 =& y_{i+1} + (y_{i+2} - y_{i+1})(t-1) \\
 *   y(x) =& s_1 + (s_2 - s_1) 0.5 t
 * \f}
 * where \f$h\f$ is the @ref step – interval between two consecutive independent variable values.
 * The function values themselves (\f$y_i\f$) are saved in the @ref params array.
 * Interpolation between the last two knots is not well-defined.
 *
 * The derivative estimation provided by @ref Diff() and @ref Diff0() is based on the formula
 * \f{equation}{
 *     y'(x)=\frac{1}{h}\left[\left(y_{i} - 2 y_{i+1} +y_{i+2} \right) t + 2 y_{i+1} - 1.5 y_{i} - 0.5 y_{i+2}\right]
 * \f}
 * which can be derived either by differentiating the interpolating formula 
 * or from the Taylor expansion of \f$f(x)\f$ around the desired independent variable value.
 *
 * The interpolator can be defined with `double` or `long double` (or even `float`) precision.
 * Error scales as \f$\mathcal{0}(h^3)\f$.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */
template <typename T>
class Linear3PInterpolator : public AbstractInterpolator<T>
{
private:
    /**
     * @brief Default constructor of Linear3PInterpolator object.
     *
     * Made private to disable its use.
     */
    Linear3PInterpolator(){};

public:
    /**
     * @brief Destroy the Linear3PInterpolator object
     *
     * Default and only destructor. Call to parent destructor should be enough.
     */
    ~Linear3PInterpolator(){};
    /**
     * @brief Construct a new Linear3PInterpolator object
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
    Linear3PInterpolator(T *knotvalues, int N, T xMax, T xMin = 0.0L, T *derivatives = nullptr)
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
     * @brief Construct a new Linear3PInterpolator object
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
    Linear3PInterpolator(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
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
     * @brief Construct a new Linear3PInterpolator object.
     *
     * Copy constructor for Linear3PInterpolator.
     *
     * @param C Original interpolator that should be copied.
     */
    Linear3PInterpolator(const Linear3PInterpolator &C)
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
     * Assign the value of @p C to the left-hand-side instance of Linear3PInterpolator.
     *
     * @note This is not a constructor!
     *
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    Linear3PInterpolator &operator=(const Linear3PInterpolator &C)
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
     * Three-point linear interpolation is used:
     * \f{aligned}{
     *   i =& \lfloor (x - x_0) / h \rfloor \\
     *   t =& \frac{x - x_i}{h} \\
     *   s_1 =& y_i + (y_{i+1} - y_i)t \\
     *   s_2 =& y_{i+1} + (y_{i+2} - y_{i+1})(t-1) \\
     *   y(x) =& s_1 + (s_2 - s_1) 0.5 t
     * \f}
     * where \f$h\f$ is the interval size between two independent variable values.
     *
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    T Eval(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue - this->step); // last interval cannot be used
        int i = (int)floor((xValue - this->minValue) / this->step);
        T t = (xValue - (this->minValue + i * this->step)) / this->step;
        T s1 = this->params[i] + (this->params[i + 1] - this->params[i]) * t;
        T s2 = this->params[i + 1] + (this->params[i + 2] - this->params[i + 1]) * (t - 1);

        return s1 + (s2 - s1) * 0.5 * t;
    };
    /**
     * @brief Evaluate function value at @p xValue.
     *
     * Evaluate the (approximation of) function value at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0 (\f$x_0 = 0.0\f$).
     * Three-point linear interpolation is used:
     * \f{aligned}{
     *   i =& \lfloor x / h \rfloor \\
     *   t =& \frac{x - ih}{h} \\
     *   s_1 =& y_i + (y_{i+1} - y_i)t \\
     *   s_2 =& y_{i+1} + (y_{i+2} - y_{i+1})(t-1) \\
     *   y(x) =& s_1 + (s_2 - s_1) 0.5 t
     * \f}
     * where \f$h\f$ is the interval size between two independent variable values.
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
        assert(xValue <= this->maxValue - this->step); // last interval cannot be used
        assert(this->minValue == 0.0);
        int i = (int)floor(xValue / this->step);
        T t = (xValue - i * this->step) / this->step;
        T s1 = this->params[i] + (this->params[i + 1] - this->params[i]) * t;
        T s2 = this->params[i + 1] + (this->params[i + 2] - this->params[i + 1]) * (t - 1);

        return s1 + (s2 - s1) * 0.5 * t;
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * The following formula to approximate the derivative is used:
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=\frac{1}{h}\left[\left(f(x_{i}) - 2 f(x_{i+1}) +f(x_{i+2}) \right) t + 2 f(x_{i+1}) - 1.5 f(x_{i}) - 0.5 f(x_{i+2})\right]
     * \f}
     * where \f$h\f$ is the @ref step and \f$t=x-x_i\f$.
     * The formula can be derived by differentiating the interpolation formula
     * or from the Taylor expansion of of \f$f(x)\f$ around \f$x=\f$ @p xValue.
     *
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    T Diff(T xValue) const override
    {
        assert(xValue >= this->minValue);
        assert(xValue <= this->maxValue - this->step);
        int i = (int)floor((xValue - this->minValue) / this->step);
        T t = (xValue - i * this->step - this->minValue) / this->step;
        return ((t - 1.5L) * this->params[i] + (2.0L - 2.0L * t) * this->params[i + 1] + (t - 0.5L) * this->params[i + 2]) / this->step;
    };
    /**
     * @brief Evaluate function derivative at @p xValue.
     *
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0 (\f$x_0 = 0.0\f$).
     * The following formula to approximate the derivative is used:
     * for \f$x\in[x_i, x_{i+1}]\f$
     * \f{equation}{
     *     f'(x)=\frac{1}{h}\left[\left(f(x_{i}) - 2 f(x_{i+1}) +f(x_{i+2}) \right) t + 2 f(x_{i+1}) - 1.5 f(x_{i}) - 0.5 f(x_{i+2})\right]
     * \f}
     * where \f$h\f$ is the @ref step and \f$t=x-x_i\f$.
     * The formula can be derived by differentiating the interpolation formula
     * or from the Taylor expansion of of \f$f(x)\f$ around \f$x=\f$ @p xValue.
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
        assert(xValue <= this->maxValue - this->step);
        int i = (int)floor(xValue / this->step);
        T t = (xValue - i * this->step) / this->step;
        return ((t - 1.5L) * this->params[i] + (2.0L - 2.0L * t) * this->params[i + 1] + (t - 0.5L) * this->params[i + 2]) / this->step;
    };
    /**
     * @brief Extend the range to zero ( @ref minValue = 0) to enable @ref Eval0().
     *
     * Extend the interpolation range to zero (set @ref minValue = 0.0) to enable Eval0().
     * The extension is possible only if @ref minValue is divisible by @ref step without remainder.
     * 3-point linear extrapolation from the first three points is used.
     *
     * @return Number of added knots (or equivalent number) or -1 if the extension cannot be done.
     */
    int ExtendToZero() override
    {
        int i;
        T *old_params = this->params;
        int no_added = (int)roundl(this->minValue / this->step);
        T t, s2, s1;
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
            // 3-point extrapolation
            t = (T)(i - no_added);
            s1 = old_params[0] + (old_params[1] - old_params[0]) * t;
            s2 = old_params[1] + (old_params[2] - old_params[1]) * (t - 1.0L);
            this->params[i] = s1 + (s2 - s1) * 0.5L * t;
            // linear extrapolation (above is better)
            // this->params[i] = old_params[0] - (old_params[1] - old_params[0]) * (no_added - i);
        }
        this->minValue = 0.0;
        delete[] old_params;

        return no_added;
    };
    void PrintInfo(std::ofstream &stream, int indentation = 0) const override
    {
        std::string indent((size_t) indentation, ' ');
        stream << indent << "- interpolation type: linear three-point interpolation (DL_POLY type, see `Linear3PInterpolator.hpp`)\n";
        stream << indent << "- usable interval: [" << this->minValue << ", " << this->maxValue << "]\n";
        stream << indent << "- resolution: " << (int) round(1.0/this->step) << " points per unit in argument (squared distance)\n";
    };
};

#endif
