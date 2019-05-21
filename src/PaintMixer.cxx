/**
 * @file PaintCoeff.cxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */

#include "PaintMixer.hxx"

#include "ColorExtraction.hxx"

#include <iomanip>
#include <stdexcept>
#include <thread>

#include <ceres/ceres.h>

namespace ceres
{
template <class T>
inline T coth(const T& x) // hyperbolic cotangent
{
  if (x == T(0))
    {
      return std::numeric_limits<T>::infinity();
    }
  else
    {
      return (ceres::exp(x) + ceres::exp(-x)) /
             (ceres::exp(x) - ceres::exp(-x));
    }
}

} // namespace ceres

namespace MixSolver
{

constexpr auto KStride = 4U;
constexpr auto Wsum    = 1.0;
constexpr auto Wsparse = 0.1;

/**
 * @brief Cost function minimizing the difference of a mixed paint from base
 * pigments and target paint.
 */
struct CostFunction_MixPaint
{
  CostFunction_MixPaint(const PaintMixer::Palette&    palette,
                        const PaintMixer::PaintCoeff& target)
    : palette(palette), target(target)
  {
  }

  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const
  {
    T Kr = T(0);
    T Kg = T(0);
    T Kb = T(0);
    T Sr = T(0);
    T Sg = T(0);
    T Sb = T(0);

    for (size_t i = 0; i < palette.size(); i++)
      {
        T c = parameters[0][i];
        Kr += c * palette[i].K[0U];
        Kg += c * palette[i].K[1U];
        Kb += c * palette[i].K[2U];
        Sr += c * palette[i].S[0U];
        Sg += c * palette[i].S[1U];
        Sb += c * palette[i].S[2U];
      }

    // lower bounds take care of negative weights.
    residuals[0] = Kr - target.K[0U];
    residuals[1] = Kg - target.K[1U];
    residuals[2] = Kb - target.K[2U];
    residuals[3] = Sr - target.S[0U];
    residuals[4] = Sg - target.S[1U];
    residuals[5] = Sb - target.S[2U];

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ::ceres::CostFunction* Create(const PaintMixer::Palette&    palette,
                                       const PaintMixer::PaintCoeff& target)
  {
    auto* c =
      (new ::ceres::DynamicAutoDiffCostFunction<CostFunction_MixPaint, KStride>(
        new CostFunction_MixPaint(palette, target)));
    c->SetNumResiduals(6);
    c->AddParameterBlock(palette.size());
    return c;
  }

  PaintMixer::Palette    palette;
  PaintMixer::PaintCoeff target;
};

/**
 * @brief Cost function for penalizing result vector whose sum is not close
 * to 1.
 */
struct CostFunction_E_sum
{
  CostFunction_E_sum(size_t n) : n(n) {}

  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const
  {

    T c0 = parameters[0][0];

    for (size_t i = 1; i < n; i++)
      {
        c0 += parameters[0][i];
      }

    residuals[0] = T(1.) - c0;

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ::ceres::CostFunction* Create(size_t n)
  {
    ::ceres::DynamicAutoDiffCostFunction<CostFunction_E_sum, KStride>* c =
      (new ::ceres::DynamicAutoDiffCostFunction<CostFunction_E_sum, KStride>(
        new CostFunction_E_sum(n)));
    c->SetNumResiduals(1);
    c->AddParameterBlock(n);
    return c;
  }

  size_t n;
};

/**
 * @brief Cost function for penalizing dense results.
 */
struct CostFunction_E_sparse
{
  CostFunction_E_sparse(size_t n) : n(n) {}

  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const
  {
    //        for (size_t i = 0; i < n; i++) {
    //            residuals[i] = T(1) - parameters[0][i];
    //        }
    //        return true;

    T l1(0);
    T l2(0);

    // https://math.stackexchange.com/questions/101200/sparseness-of-a-vector
    for (size_t i = 0; i < n; i++)
      {
        l1 += parameters[0][i];
        l2 += ceres::pow(parameters[0][i], T(2));
      }
    l2 = ceres::sqrt(l2);

    T s = ceres::sqrt(T(n));

    residuals[0] = T(1) - ((s - (l1 / l2)) / (s - T(1)));

    //        std::cout << "Density: " <<std::setprecision(3) <<residuals[0]  <<
    //        std::endl;

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ::ceres::CostFunction* Create(size_t n)
  {
    ::ceres::DynamicAutoDiffCostFunction<CostFunction_E_sparse, KStride>* c =
      (new ::ceres::DynamicAutoDiffCostFunction<CostFunction_E_sparse, KStride>(
        new CostFunction_E_sparse(n)));
    c->SetNumResiduals(1);
    c->AddParameterBlock(n);
    return c;
  }

  size_t n;
};

} // namespace MixSolver

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
 * @brief Get the mixture recipe for the underlying palette to mix a target
 * paint.
 *
 * @param paint the target paint.
 * @return std::vector<CoeffPrecision>
 */
std::vector<CoeffPrecision>
PaintMixer::getWeightsForMixingTargetPaint(const PaintCoeff& paint) const
{
  const auto k = mBasePalette.size();

  std::vector<CoeffPrecision> weights(k);

  for (auto i = 0U; i < k; i++)
    {
      weights[i] = 1.0 / static_cast<float64_t>(k);
    }

  ceres::Problem problem;

  ::ceres::CostFunction* dataCostFunction =
    MixSolver::CostFunction_MixPaint::Create(mBasePalette, paint);
  problem.AddResidualBlock(dataCostFunction, nullptr, &(weights.front()));

  ::ceres::CostFunction* sumCostFunction =
    MixSolver::CostFunction_E_sum::Create(k);
  problem.AddResidualBlock(
    sumCostFunction,
    new ceres::ScaledLoss(nullptr, MixSolver::Wsum, ceres::TAKE_OWNERSHIP),
    weights.data());

  ::ceres::CostFunction* sumSparseFunction =
    MixSolver::CostFunction_E_sparse::Create(k);
  problem.AddResidualBlock(
    sumSparseFunction,
    new ceres::ScaledLoss(nullptr, MixSolver::Wsparse, ceres::TAKE_OWNERSHIP),
    weights.data());

  problem.AddResidualBlock(dataCostFunction, nullptr, &(weights.front()));
  for (auto i = 0U; i < k; ++i)
    {
      problem.SetParameterLowerBound(&(weights.front()), i, 0.0);
      problem.SetParameterUpperBound(&(weights.front()), i, 1.0);
    }

  ::ceres::Solver::Options options;
  //    options.minimizer_progress_to_stdout = true;
  const auto nThreads =
    static_cast<int32_t>(std::thread::hardware_concurrency());
  options.num_threads               = nThreads;
  options.num_linear_solver_threads = nThreads;
  options.max_num_iterations        = 100;
  options.function_tolerance        = 1e-6;

  ::ceres::Solver::Summary summary;
  ::ceres::Solve(options, &problem, &summary);
  //    LOG(INFO) << summary.BriefReport() << "\n";

  float64_t         wSum = 0.0;
  std::stringstream stream;
  stream << "weights: ";
  for (auto i = 0U; i < weights.size(); i++)
    {
      stream << std::setprecision(3) << weights[i] << "\t";

      wSum += weights[i];
    }
  stream << "| sum: " << std::setprecision(3) << wSum << std::endl;
  LOG(INFO) << stream.str();
  if (!fuzzyEqual(wSum, 1.0))
    {
      // normalize to sum one
      float64_t norm = 1.0 / wSum;
      for (size_t i = 0; i < weights.size(); i++)
        {
          weights[i] *= norm;
        }
    }
  return weights;
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
