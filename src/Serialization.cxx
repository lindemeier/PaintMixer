/**
 * @file Serialization.cxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */

#include "PaintMixer/Serialization.hxx"
#include "PaintMixer/KubelkaMunk.hxx"

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <ColorConverter/ColorConverter.hxx>

namespace PaintMixer
{
void LoadPalette(std::istream& stream, Palette& palette)
{
  cereal::JSONInputArchive ar(stream);
  ar(cereal::make_nvp("palette", palette));
}

void SavePalette(std::ostream& stream, const Palette& palette)
{
  cereal::JSONOutputArchive ar(stream);
  ar(cereal::make_nvp("palette", palette));
}

std::string ExtractFiletype(const std::string& filename)
{
  std::string res(filename);
  size_t      ls = res.find_last_of(".");
  res            = res.substr(ls + 1, res.size() - ls - 1);
  return res;
}

/**
 * @brief Load a picture from a file.
 *
 * @param filenameOriginal
 * @param image
 */
void LoadImage(const std::string& filenameOriginal, cv::Mat_<vec3f>& image)
{
  std::string filename = filenameOriginal;
  std::replace(filename.begin(), filename.end(), '\\', '/');

  cv::Mat cv_mat =
    cv::imread(filename, cv::IMREAD_ANYDEPTH | cv::IMREAD_COLOR);

  // if not loaded succesfully
  if (!cv_mat.data)
    {
      throw std::ios_base::failure(filenameOriginal);
    }

  if (cv_mat.channels() == 1)
    {
      cv::Mat in[] = {cv_mat, cv_mat, cv_mat};
      cv::merge(in, 3, cv_mat);
    }
  else if (cv_mat.channels() == 4)
    {
      cv::cvtColor(cv_mat, cv_mat, cv::COLOR_BGRA2BGR);
    }

  // data scale
  float32_t scale = 1.0f;
  if (cv_mat.depth() == CV_16U)
    scale = 1.0f / (0xffff);
  else if (cv_mat.depth() == CV_32F)
    scale = 1.0f;
  else if (cv_mat.depth() == CV_8U)
    scale = 1.0f / (0xff);
  else if (cv_mat.depth() == CV_64F)
    scale = 1.0f / (0xffffffff);

  // convert to right type
  cv_mat.convertTo(cv_mat, CV_32FC3, scale);

  // OpenCV has BGR
  cv::cvtColor(cv_mat, image, cv::COLOR_BGR2RGB);
}

/**
 * @brief Save a picture to file.
 *
 * @param filenameOriginal
 * @param output
 *
 * @return true
 * @return false
 */
bool SaveImage(const std::string&     filenameOriginal,
               const cv::Mat_<vec3f>& output)
{
  std::string filename = filenameOriginal;
  std::replace(filename.begin(), filename.end(), '\\', '/');

  std::string filetype = ExtractFiletype(filename);

  cv::Mat m;
  if (filetype == "png")
    {
      const auto scale = static_cast<float32_t>(0xffff);
      output.convertTo(m, CV_MAKETYPE(CV_16U, 3), scale);
    }
  else
    {
      const auto scale = static_cast<float32_t>(0xff);
      output.convertTo(m, CV_MAKETYPE(CV_8U, 3), scale);
    }
  cv::cvtColor(m, m, cv::COLOR_RGB2BGR);
  std::vector<int32_t> params;
  params.push_back(cv::IMWRITE_JPEG_QUALITY);
  params.push_back(100);
  params.push_back(cv::IMWRITE_PNG_COMPRESSION);
  params.push_back(0);
  return cv::imwrite(filename, m, params);
}

/**
 * @brief Visualize the palette colors by applying them onto black and white
 * backgrounds.
 *
 * @param palette
 * @param appliedThickness
 *
 * @return cv::Mat_<vec3f>
 */
cv::Mat_<vec3f> VisualizePalette(const Palette&  palette,
                                 const float64_t appliedThickness)
{
  cv::Mat_<vec3f> paletteImage(200, palette.size() * 100);
  using Converter = color::ColorConverter<float64_t, 3UL, vec>;
  Converter converter;

  vec3d black = {0., 0., 0.};
  vec3d white = {1., 1., 1.};

  for (size_t c = 0; c < palette.size(); c++)
    {
      const auto& paint = palette[c];

      vec3d colorOnBlackRGB =
        ComputeReflectance(paint, black, appliedThickness);
      vec3d colorOnWhiteRGB =
        ComputeReflectance(paint, white, appliedThickness);

      vec3d colorOnBlack, colorOnWhite;
      converter.rgb2srgb(colorOnBlackRGB, colorOnBlack);
      converter.rgb2srgb(colorOnWhiteRGB, colorOnWhite);

      const int32_t h = paletteImage.rows;

      for (uint32_t x = c * (paletteImage.cols / palette.size());
           x < (c + 1) * (paletteImage.cols / palette.size()); x++)
        {
          int32_t st = 0;

          for (int32_t y = st; y < 0.1 * h; y++)
            {
              paletteImage(y, x)[0] = 1.f;
              paletteImage(y, x)[1] = 1.f;
              paletteImage(y, x)[2] = 1.f;
            }
          st = 0.1 * h;

          for (int32_t y = st; y < 0.4 * h; y++)
            {
              paletteImage(y, x)[0] = colorOnWhite[0];
              paletteImage(y, x)[1] = colorOnWhite[1];
              paletteImage(y, x)[2] = colorOnWhite[2];
            }
          st = 0.4 * h;

          for (int32_t y = st; y < 0.5 * h; y++)
            {
              paletteImage(y, x)[0] = 1.f;
              paletteImage(y, x)[1] = 1.f;
              paletteImage(y, x)[2] = 1.f;
            }
          st = 0.5 * h;

          for (int32_t y = st; y < 0.6 * h; y++)
            {
              paletteImage(y, x)[0] = 0.f;
              paletteImage(y, x)[1] = 0.f;
              paletteImage(y, x)[2] = 0.f;
            }
          st = 0.6 * h;

          for (int32_t y = st; y < 0.9 * h; y++)
            {
              paletteImage(y, x)[0] = colorOnBlack[0];
              paletteImage(y, x)[1] = colorOnBlack[1];
              paletteImage(y, x)[2] = colorOnBlack[2];
            }
          st = 0.9 * h;

          for (int32_t y = st; y < h; y++)
            {
              paletteImage(y, x)[0] = 0.f;
              paletteImage(y, x)[1] = 0.f;
              paletteImage(y, x)[2] = 0.f;
            }
        }
    }
  return paletteImage;
}

} // namespace PaintMixer
