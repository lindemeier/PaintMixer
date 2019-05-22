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
#include <iostream>

#include "Types.hxx"

namespace PaintMixer
{

constexpr auto CoeffSamplesCount = 3U;
using CoeffPrecision             = float64_t;

/**
 * @brief Represents Kubelka-Munk coefficients.
 *
 */
struct PaintCoeff
{
  using VecType = std::array<CoeffPrecision, CoeffSamplesCount>;

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

std::ostream& operator<<(std::ostream& output, const PaintMixer::PaintCoeff& v);

#endif // PAINT_MIXER_PAINT_COEFF_H
