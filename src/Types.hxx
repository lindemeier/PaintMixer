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

#include <Eigen/Core>
#include <opencv2/core.hpp>

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

namespace cv
{ // opencv access data traits
template <typename T, int32_t R, int32_t C>
class DataType<Eigen::Matrix<T, R, C>>
{
public:
  typedef Eigen::Matrix<T, R, C> value_type;
  typedef Eigen::Matrix<T, R, C> work_type;
  typedef T                      channel_type;
  typedef value_type             vec_type;
  enum
  {
    generic_type = 0,
    depth        = DataDepth<channel_type>::value,
    channels     = R * C,
    fmt          = ((channels - 1) << 8) + DataDepth<channel_type>::fmt,
    type         = CV_MAKETYPE(depth, channels)
  };
};

template <typename T, int32_t R>
class DataType<std::array<T, R>>
{
public:
  typedef std::array<T, R> value_type;
  typedef std::array<T, R> work_type;
  typedef T                channel_type;
  typedef value_type       vec_type;
  enum
  {
    generic_type = 0,
    depth        = DataDepth<channel_type>::value,
    channels     = R,
    fmt          = ((channels - 1) << 8) + DataDepth<channel_type>::fmt,
    type         = CV_MAKETYPE(depth, channels)
  };
};

} // namespace cv

#endif // PAINT_MIXER_TYPES_H