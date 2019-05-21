/**
 * @file PaintCoeff.cxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */

#include "PaintMixer.hxx"

#include "ColorExtraction.hxx"

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
 * @param sRGBPicture The input image. This should be in linear rgb.
 * @param count The number of paints in the resulting palette.
 * @return Palette the palette
 */
Palette PaintMixer::mixFromInputPicture(const cv::Mat_<vec3f>& sRGBPicture,
                                        uint32_t               count) const
{
  // extract rgb palette from the image
  std::vector<vec3f> colors;
  ExtractColorPaletteAharoni(sRGBPicture, colors, count);

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
  if (!fuzzyEqual(wSum, static_cast<float64_t>(1.0)))
    {
      // normalize to sum one
      norm = 1. / wSum;
    }

  PaintCoeff p;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      p.K[i] = 0.0;
      p.S[i] = 0.0;
    }

  for (auto l = 0U; l < mBasePalette.size(); l++)
    {
      for (auto i = 0U; i < CoeffSamplesCount; i++)
        {
          p.K[i] += norm * weights[l] * mBasePalette[l].K[i];
          p.S[i] += norm * weights[l] * mBasePalette[l].S[i];
        }
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
