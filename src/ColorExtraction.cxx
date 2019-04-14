/**
 * @file ColorExtraction.hxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-14
 */
#include "ColorExtraction.hxx"

#include <opencv2/imgproc.hpp>

namespace PaintMixer
{
/**
 * Pigment-Based Recoloring of Watercolor Paintings -
 * Elad Aharoni-Mack, Yakov Shambik and Dani Lischinski -
 * Expressive 2017
 * @param sRGB input picture in srgb space
 * @param colors output vector with RGB colors
 * @param k number of colors
 */
void ExtractColorPaletteAharoni(
  const cv::Mat_<color::ColorConverter<float>::vec3>& sRGB,
  std::vector<color::ColorConverter<float>::vec3>& linearRGB_colors, int32_t k)
{
  using vec3  = color::ColorConverter<float>::vec3;
  using vec2f = Eigen::Matrix<float, 2U, 1U>;

  // https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
  const auto DistanceLinePoint = [](const vec2f& v, const vec2f& w,
                                    const vec2f& p) {
    // Return minimum distance between line segment vw and point p
    const auto l2 = (v - w).squaredNorm(); // i.e. |w-v|^2 -  avoid a sqrt
    if (std::fabs(l2) < std::numeric_limits<float>::epsilon())
      {
        return (p - v).norm(); // v == w case
      }

    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line.
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    // We clamp t from [0,1] to handle points outside the segment vw.
    const auto t =
      std::max<float>(0., std::min<float>(1., (p - v).dot(w - v) / l2));
    const vec2f projection = v + t * (w - v); // Projection falls on the segment

    return (p - projection).norm();
  };

  // convert to lab space
  std::vector<vec3>            Lab(sRGB.total());
  color::ColorConverter<float> converter;

  for (auto i = 0U; i < Lab.size(); i++)
    {
      converter.convert(
        sRGB(i), Lab[i],
        color::ColorConverter<float>::Conversion::srgb_2_CIELab);
    }

  auto       L_max  = 0.f;
  auto       L_min  = 100.f;
  const auto L_down = 0.f;
  const auto L_up   = 0.f;
  vec3       c_max_L;
  vec3       c_min_L;
  for (const vec3& e : Lab)
    {
      if (L_max < e[0])
        {
          L_max   = e[0];
          c_max_L = e;
        }
      if (L_min > e[0])
        {
          L_min   = e[0];
          c_min_L = e;
        }
    }

  // list of colors in ab plane
  std::vector<vec2f>   inputPoints;
  std::vector<vec3>    inputColors;
  std::vector<int32_t> indices;
  inputPoints.reserve(Lab.size());
  inputColors.reserve(Lab.size());
  indices.reserve(Lab.size());
  int32_t index = 0;
  for (const vec3& e : Lab)
    {
      // To avoid noise, we first discard the brightest and the darkest
      // colors, and compute the convex hull of the remaining colors in the
      // image.
      if (e[0] > (L_min + L_down) && e[0] < (L_max - L_up))
        {
          inputPoints.push_back({e[1], e[2]});
          inputColors.push_back(e);
          indices.push_back(index++);
        }
    }

  indices.clear();
  cv::convexHull(inputPoints, indices, true, false);

  linearRGB_colors.reserve(k);
  k -= 2;
  vec3 a, b;
  converter.convert(c_min_L, a,
                    color::ColorConverter<float>::Conversion::CIELab_2_rgb);
  linearRGB_colors.push_back(a);
  converter.convert(c_max_L, b,
                    color::ColorConverter<float>::Conversion::CIELab_2_rgb);
  linearRGB_colors.push_back(b);

  // trivial cases
  if (static_cast<int32_t>(indices.size()) <= k ||
      k > static_cast<int32_t>(indices.size()))
    {
      for (auto i : indices)
        {
          vec3 a;
          converter.convert(
            inputColors[i], a,
            color::ColorConverter<float>::Conversion::CIELab_2_rgb);
          linearRGB_colors.push_back({a[0], a[1], a[2]});
        }
      return;
    }
  /*
   * We then greedily iterate in order to simplify the convex hull
   * polygon by pruning its vertices until we are left with k vertices,
   * where k is the desired palette size. Pruning is done similarly to
   * the Douglas-Peucker algorithm [1973], where at each iteration
   * we remove the vertex whose distance from the line connecting
   * its neighbors is the smallest.
   * */
  while (static_cast<int32_t>(indices.size()) > k)
    {
      // first point
      vec2f p_1 = inputPoints[indices.back()];
      vec2f p0  = inputPoints[indices.front()];
      vec2f p1  = inputPoints[indices[1]];

      size_t removeIndex = 0;
      float  md          = DistanceLinePoint(p_1, p1, p0);

      // inner points
      for (size_t l = 1; l < indices.size() - 1; l++)
        {
          const int32_t i_1 = indices[l - 1];
          const int32_t i0  = indices[l];
          const int32_t i1  = indices[l + 1];

          p_1 = inputPoints[i_1];
          p0  = inputPoints[i0];
          p1  = inputPoints[i1];

          float d = DistanceLinePoint(p_1, p1, p0);
          if (d < md)
            {
              md          = d;
              removeIndex = l;
            }
        }

      // last point
      p_1     = inputPoints[*(indices.cend() - 2)];
      p0      = inputPoints[indices.back()];
      p1      = inputPoints[indices.front()];
      float d = DistanceLinePoint(p_1, p1, p0);
      if (d < md)
        {
          removeIndex = indices.size() - 1;
        }

      indices.erase(indices.begin() + removeIndex);
    }

  for (auto i : indices)
    {
      vec3 a;
      converter.convert(inputColors[i], a,
                        color::ColorConverter<float>::Conversion::CIELab_2_rgb);
      linearRGB_colors.push_back({a[0], a[1], a[2]});
    }
}
} // namespace PaintMixer