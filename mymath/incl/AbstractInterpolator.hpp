/*
 *  Interface for classes providing function evaluation through interpolation and splines.
 *  Originally for Simul++ simulation package
 *
 *  Author: JJ
 *  Date: March 2023
 *
 */

#ifndef ABSTRACTINTERPOLATORHEADER
#define ABSTRACTINTERPOLATORHEADER

/**
 * @defgroup interp Function interpolators
 * @ingroup math
 * @brief Set of classes that provide complicated function evaluation by interpolation and/or splines.
 *
 * Calculation of a complicated function can cost a lot of CPU time.
 * Classes in this group provide faster evaluation of functions using their approximation by various types of interpolation.
 *
 *
 * @class AbstractInterpolator
 * @ingroup interp
 * @brief Interface for classes providing function evaluation through interpolation and splines.
 *
 * This class provides interface for creation and use of interpolators used to evaluate complicated
 * functions using interpolation and/or splines.
 * Once initialized, the interpolator should be able to compute function values by @ref Eval() and @ref Eval0()
 * much faster than if they were calculated using direct formulae.
 * Moreover, the derivative can be estimated by @ref Diff() or @ref Diff0().
 * The interface can be defined with `double` or `long double` precision by providing template @p T.
 *
 * @tparam T `double` or `long double` (or `float`). Determines the precision of all calculations inside.
 */
template <typename T>
class AbstractInterpolator
{
protected:
    /**
     * @brief Array of parameters needed to calculate the interpolated value
     *
     * This array is constructed once and than used for function evaluation by @ref Eval() and @ref Eval0().
     */
    T *params;
    /**
     * @brief Number of points (knots) which are interpolated.
     *
     * This number does not need to match the size of @ref params.
     * If 4 parameters are needed inside one interval to calculate the function value, the size of @ref params is 4*@ref size.
     */
    int size;
    ///@{
    /**
     * @brief Minimum and maximum allowed value of independent variable.
     *
     * The interpolator is able to compute function values only in this range.
     * Frequently, @ref minValue = 0.0 and @ref Eval0() can be used instead of @ref Eval().
     */
    T minValue;
    T maxValue;
    ///@}
    /**
     * @brief The size of the interval between two consequent knots.
     *
     * Needed to calculate the appropriate indeces in @ref params for function evaluation.
     */
    T step;
    /**
     * @brief Default constructor of Abstract Interpolator object.
     *
     * Made private to disable its use.
     */
    AbstractInterpolator(){params = nullptr;};

public:
    /**
     * @brief Destroy the Abstract Interpolator object
     *
     * Default and only destructor. Deallocates parameter array @ref params.
     */
    virtual ~AbstractInterpolator()
    {
        if (params != nullptr)
        {
            delete[] params;
        }
        params = nullptr;
    };
    /**
     * @brief Construct a new Abstract Interpolator object
     * 
     * Define function values at knots in @p knotvalues and the number of knots in @p N.
     * The interpolation range is given by @p xMin and @p xMax; the former is assumed to be 0.0 if not specified otherwise.
     * Some interpolators (splines) need values of function derivative (provided by @p derivatives) at end points or at all knots.
     * 
     * @param knotvalues Function values at equidistant knots.
     * @param N Number of knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives Function derivatives at end points or at all knots (interpolator dependent).
     */
    AbstractInterpolator(T* knotvalues, int N, T xMax, T xMin = 0.0L, T* derivatives = nullptr)
    {
        minValue = xMin;
        maxValue = xMax;
        size = N;
        params = nullptr;
    };
    /**
     * @brief Construct a new Abstract Interpolator object
     * 
     * Define function values at knots in @p knotvalues; the number of knots can be determined by @p knotvalues size.
     * The interpolation range is given by @p xMin and @p xMax; the former is assumed to be 0.0 if not specified otherwise.
     * Some interpolators (splines) need values of function derivative (provided by @p derivatives) at end points or at all knots.
     * 
     * @param knotvalues Function values at equidistant knots.
     * @param xMax Maximum independent variable value.
     * @param xMin Minimum independent variable value.
     * @param derivatives Function derivatives at end points or at all knots (interpolator dependent).
     */
    AbstractInterpolator(std::vector<T> knotvalues, T xMax, T xMin = 0.0L, std::vector<T> derivatives = std::vector<T>())
    {
        minValue = xMin;
        maxValue = xMax;
        size = knotvalues.size();
        params = nullptr;
    };
    /**
     * @brief Construct a new Abstract Interpolator object.
     * 
     * Copy constructor for AbstractInterpolator.
     * 
     * @param C Original interpolator that should be copied.
     */
    AbstractInterpolator(const AbstractInterpolator &C) 
    {
        params = nullptr; // interpolator specific size cannot be predetermined
        size = C.size;
        minValue = C.minValue;
        maxValue = C.maxValue;
        step = C.step;
    };
    /**
     * @brief Default assignment operator.
     * 
     * Assign the value of @p C to the left-hand-side instance of AbstractInterpolator.
     * 
     * @note This is not a constructor!
     * 
     * @param C The original interpolator.
     * @return Assignment copy...
     */
    virtual AbstractInterpolator &operator=(const AbstractInterpolator &C)
    {
        params = nullptr; // interpolator specific size cannot be predetermined
        size = C.size;
        minValue = C.minValue;
        maxValue = C.maxValue;
        step = C.step;
        return *this;
    }
    /**
     * @brief Evaluate function value at @p xValue.
     * 
     * Evaluate the (approximation of) function value at the point @p xValue.
     * Pure virtual to be implemented be the child class.
     * 
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    virtual T Eval(T xValue) const = 0;
     /**
     * @brief Evaluate function value at @p xValue.
     * 
     * Evaluate the (approximation of) function value at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * Pure virtual to be implemented be the child class.
     * 
     * @note No checks are provided that @ref minValue = 0.0!
     * Use to your own risk if you need to speed up the calculation a little.
     * 
     * @param xValue Independent variable value.
     * @return T Function value at @p xValue.
     */
    virtual T Eval0(T xValue) const = 0;
    /**
     * @brief Evaluate function derivative at @p xValue.
     * 
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * This version assumes that the @ref minValue = 0.0.
     * Pure virtual to be implemented be the child class.
     * 
     * @note No checks are provided that @ref minValue = 0.0!
     * Use to your own risk if you need to speed up the calculation a little.
     * 
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    virtual T Diff0(T xValue) const = 0;
    /**
     * @brief Evaluate function derivative at @p xValue.
     * 
     * Evaluate the (approximation of) function derivative at the point @p xValue.
     * Pure virtual to be implemented be the child class.
     * 
     * @param xValue Independent variable value.
     * @return T Function derivative approximation at @p xValue.
     */
    virtual T Diff(T xValue) const = 0;
    /**
     * @brief Extend the range to zero ( @ref minValue = 0) to enable @ref Eval0().
     * 
     * Extend the interpolation range to zero (set @ref minValue = 0.0) to enable Eval0().
     * Uses some heuristic to provide this **extrapolation**.
     * 
     * @return Number of added knots (or equivalent number) or -1 if the extension cannot be done.
     */
    virtual int ExtendToZero() = 0;
    /**
     * @brief Print basic information about the interpolator.
     * 
     * Used to report the minimum and maximum possible values, type of interpolation and gridsize.
     * Info is printed in a `YAML`-based format and indented by @p indentation number of whitespaces.
     * 
     * @param stream Output stream.
     * @param indentation Number of leading whitespaces to enable indentation
     */
    virtual void PrintInfo(std::ofstream& stream, int indentation = 0) const = 0;
};

#endif
