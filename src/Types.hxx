/**
 * @file Types.hxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */
#ifndef PAINT_MIXER_TYPES_H
#define PAINT_MIXER_TYPES_H

#include <cmath>
#include <stdint.h>

// TODO find a platform independent solution
using float64_t = double;
using float32_t = float;

/**
 * @brief Fuzzy comparison of floating points.
 *
 * @tparam FloatingType the floating point type.
 * @param[in] firstValue first value for comparison.
 * @param[in] secondValue second value for comparison.
 * @param[in] epsilon the fuzzyiness of the comparison.
 * @return if the given values are fuzzy equal.
 */
template <typename FloatingType>
typename std::enable_if<std::is_floating_point<FloatingType>::value, bool>::type
fuzzyEqual(
  const FloatingType firstValue, const FloatingType secondValue,
  const FloatingType epsilon = std::numeric_limits<FloatingType>::epsilon())
{
  return std::fabs(firstValue - secondValue) < epsilon;
}

#endif // PAINT_MIXER_TYPES_H