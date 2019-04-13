/**
 * @file PaintCoeff.hxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */
#ifndef PAINT_MIXER_PAINT_COEFF_H
#define PAINT_MIXER_PAINT_COEFF_H

#include <array>

#include <Eigen/Core>

#include "Types.hxx"

namespace PaintMixer
{

constexpr auto CoeffSamplesCount = 7U;
using CoeffPrecision             = float64_t;

/**
 * @brief Represents Kubelka-Munk coefficients.
 *
 */
struct PaintCoeff
{
  using VecType = Eigen::Matrix<CoeffPrecision, CoeffSamplesCount, 1U>;

  /**
   * @brief Absorption
   */
  VecType K;
  /**
   * @brief Scattering
   */
  VecType S;
};

} // namespace PaintMixer

#endif // PAINT_MIXER_PAINT_COEFF_H
