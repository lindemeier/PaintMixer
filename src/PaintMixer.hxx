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

namespace PaintMixer
{

/**
 * @brief Represents a class that offers paint mix functions.
 *
 */
class PaintMixer
{
  /**
   * @brief The underlying collection of base paints that can be used to mix
   * other Paints and Palettes.
   *
   */
  Palette mBasePalette;

public:
  PaintMixer(const Palette& basePalette);

  Palette mixFromInputPicture(
    const std::vector<std::array<float32_t, 3U>>& linearRGBPicture,
    uint32_t                                      count) const;

  PaintCoeff mixSinglePaint(const std::vector<CoeffPrecision>& weights) const;

  const Palette& getUnderlyingPalette() const;
  void           setUnderlyingPalette(const Palette& palette);

private:
  PaintMixer() = delete;
};

} // namespace PaintMixer

#endif // PAINT_MIXER_PAINT_MIXER_H