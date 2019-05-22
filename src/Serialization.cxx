/**
 * @file Serialization.cxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */

#include "PaintMixer/Serialization.hxx"

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace PaintMixer
{
void LoadPalette(std::istream& stream, Palette& palette)
{
  try
    {
      cereal::JSONInputArchive ar(stream);
      ar(cereal::make_nvp("palette", palette));
    }
  catch (cereal::Exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
}

void SavePalette(std::ostream& stream, const Palette& palette)
{
  try
    {
      cereal::JSONOutputArchive ar(stream);
      ar(cereal::make_nvp("palette", palette));
    }
  catch (cereal::Exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
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
    cv::imread(filename, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_COLOR);

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
    scale = 1.0f / (0xffff - 1);
  else if (cv_mat.depth() == CV_32F)
    scale = 1.0f;
  else if (cv_mat.depth() == CV_8U)
    scale = 1.0f / (0xff - 1);
  else if (cv_mat.depth() == CV_64F)
    scale = 1.0f / (0xffffffff - 1);

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
      const auto scale = static_cast<float32_t>(0xffff - 1);
      output.convertTo(m, CV_MAKETYPE(CV_16U, 3), scale);
    }
  else
    {
      const auto scale = static_cast<float32_t>(0xff - 1);
      output.convertTo(m, CV_MAKETYPE(CV_8U, 3), scale);
    }
  cv::cvtColor(m, m, cv::COLOR_RGB2BGR);
  std::vector<int32_t> params;
  params.push_back(CV_IMWRITE_JPEG_QUALITY);
  params.push_back(100);
  params.push_back(CV_IMWRITE_PNG_COMPRESSION);
  params.push_back(0);
  return cv::imwrite(filename, m, params);
}

} // namespace PaintMixer
