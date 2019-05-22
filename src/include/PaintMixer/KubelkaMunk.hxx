/**
 * @file KubelkaMunk.hxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-05-22
 *
 */
#ifndef PAINT_MIXER_KUBELKA_MUNK_H
#define PAINT_MIXER_KUBELKA_MUNK_H

#include "PaintMixer/PaintCoeff.hxx"
#include "PaintMixer/Types.hxx"

#include <limits>
#include <type_traits>

namespace PaintMixer
{
template <class T,
          typename std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
inline T coth(const T& x) // hyperbolic cotangent
{
  if (fuzzyEqual(x, T(0)))
    {
      return std::numeric_limits<T>::infinity();
    }
  else
    {
      return (std::exp(x) + std::exp(-x)) / (std::exp(x) - std::exp(-x));
    }
}

template <class Scalar, typename std::enable_if_t<
                          std::is_floating_point<Scalar>::value, int> = 0>
vec<Scalar, CoeffSamplesCount>
ComputeReflectance(const PaintCoeff&                     coeffs,
                   const vec<Scalar, CoeffSamplesCount>& R0, const Scalar d)
{
  using vec_type = vec<Scalar, CoeffSamplesCount>;

  if (d <= 0.0)
    {
      return R0;
    }
  vec_type K_S;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      K_S[i] = coeffs.K[i] / coeffs.S[i];
    }

  vec_type a;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      a[i] = Scalar(1.0) + K_S[i];
    }

  vec_type asq;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      asq[i] = a[i] * a[i];
    }

  vec_type b;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      b[i] = std::sqrt(asq[i] - Scalar(1.0));
    }

  vec_type bSh;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      bSh[i] = b[i] * coeffs.S[i] * d;
    }

  vec_type bcothbSh;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      bcothbSh[i] = b[i] * coth(bSh[i]);
    }
  vec_type R;
  for (auto i = 0U; i < CoeffSamplesCount; i++)
    {
      R[i] = (Scalar(1.0) - R0[i] * (a[i] - bcothbSh[i])) /
             (a[i] - R0[i] + bcothbSh[i]);
    }

  return R;
}

} // namespace PaintMixer

#endif // PAINT_MIXER_KUBELKA_MUNK_H