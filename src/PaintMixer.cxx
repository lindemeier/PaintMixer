/**
 * @file PaintCoeff.cxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */

#include "PaintMixer.hxx"

#include <stdexcept>

namespace PaintMixer
{
/**
 * @brief Construct a new Paint Mixer::Paint Mixer object
 *
 * @param basePalette the underlying base palette.
 */
PaintMixer::PaintMixer(const Palette& basePalette) : mBasePalette(basePalette)
{
}

/**
 * @brief Mix a palette from an input RGB image. The image is analyzed using a
 * color clustering algorithm. The best count fits are used to mix a fitting
 * palette from the base palette.
 *
 * @param linearRGBPicture The input image. This should be in linear rgb.
 * @param count The number of paints in the resulting palette.
 * @return Palette the palette
 */
Palette PaintMixer::mixFromInputPicture(
  const std::vector<std::array<float32_t, 3U>>& linearRGBPicture,
  uint32_t                                      count) const
{
  // TODO
  return Palette();
}

/**
 * @brief Mix a single paint from the base palette according to the given
 * weights. (Weighted linear combination)
 *
 * @param weights The weights used for mixing.
 * @return PaintCoeff
 */
PaintCoeff
PaintMixer::mixSinglePaint(const std::vector<CoeffPrecision>& weights) const
{
  if (weights.size() != mBasePalette.size())
    {
      throw std::invalid_argument(
        "Palette size does not match underlying size.");
    }
  auto norm = static_cast<float64_t>(1.0);
  auto wSum = static_cast<float64_t>(0.0);
  for (size_t i = 0; i < weights.size(); i++)
    {
      wSum += weights[i];
    }
  if (!fuzzyEqual(wSum, static_cast<float64_t>(1.0), 0.000001))
    {
      // normalize to sum one
      norm = 1. / wSum;
    }

  PaintCoeff p = {PaintCoeff::VecType::Zero(), PaintCoeff::VecType::Zero()};

  for (auto i = 0U; i < mBasePalette.size(); i++)
    {
      p.K += norm * weights[i] * mBasePalette[i].K;
      p.S += norm * weights[i] * mBasePalette[i].S;
    }
  return p;
}

/**
 * @brief Get the Palette object
 *
 * @return const Palette&
 */
const Palette& PaintMixer::getUnderlyingPalette() const { return mBasePalette; }

/**
 * @brief Set the Palette object
 *
 * @param palette
 */
void PaintMixer::setUnderlyingPalette(const Palette& palette)
{
  mBasePalette = palette;
}

} // namespace PaintMixer
