/**
 * @file PaintMixer.hxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */
#ifndef PAINT_MIXER_PAINT_MIXER_H
#define PAINT_MIXER_PAINT_MIXER_H

#include "Palette.hxx"

#include <opencv2/core.hpp>

namespace PaintMixer
{

/**
 * @brief Represents a class that offers paint mix functions.
 *
 */
class PaintMixer
{

public:
  PaintMixer(const Palette& basePalette);

  Palette mixFromInputPicture(const cv::Mat_<vec3f>& sRGBPicture,
                              uint32_t               count) const;

  PaintCoeff mixSinglePaint(const std::vector<CoeffPrecision>& weights) const;

  std::vector<CoeffPrecision>
  getWeightsForMixingTargetPaint(const PaintCoeff& paint) const;

  std::vector<CoeffPrecision>
  getMixtureWeightsForReflectance(const vec3d& targetReflectance,
                                  const vec3d& backgroundReflectance,
                                  float64_t&   layerThickness) const;

  const Palette& getUnderlyingPalette() const;
  void           setUnderlyingPalette(const Palette& palette);

private:
  PaintMixer() = delete;

  /**
   * @brief The underlying collection of base paints that can be used to mix
   * other Paints and Palettes.
   *
   */
  Palette mBasePalette;
};

} // namespace PaintMixer

#endif // PAINT_MIXER_PAINT_MIXER_H