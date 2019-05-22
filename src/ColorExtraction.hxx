/**
 * @file ColorExtraction.hxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-14
 */
#ifndef PAINT_MIXER_COLOR_EXTRACTION_H
#define PAINT_MIXER_COLOR_EXTRACTION_H

#include "PaintMixer/Types.hxx"

namespace PaintMixer
{

void ExtractColorPaletteAharoni(const cv::Mat_<vec3f>& sRGB,
                                std::vector<vec3f>&    linearRGB_colors,
                                int32_t                k);
} // namespace PaintMixer
#endif // PAINT_MIXER_COLOR_EXTRACTION_H