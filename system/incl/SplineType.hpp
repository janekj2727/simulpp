#ifndef SPLINETYPEENUM
#define SPLINETYPEENUM

/**
 * @brief Enum for type of splines.
 * 
 * Used during initialization of splines to be readable.
 */
enum SplineType
{
  hyperbolic,
  quadratic,
  natural,
  hermite,
  linear,
  linear3,
  linear4
};
#endif