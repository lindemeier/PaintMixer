/**
 * @file PaintCoeff.cxx
 * @author Thomas Lindemeier
 * @brief
 * @date 2019-04-13
 *
 */

#include "PaintCoeff.hxx"

namespace PaintMixer
{

} // namespace PaintMixer

std::ostream& operator<<(std::ostream& output, const PaintMixer::PaintCoeff& v)
{
  output << "K(";
  size_t i = 0;
  for (; i < PaintMixer::CoeffSamplesCount - 1; i++)
    {
      output << v.K[i] << ", ";
    }
  output << v.K[i] << ")";
  output << "\tS(";
  i = 0;
  for (; i < PaintMixer::CoeffSamplesCount - 1; i++)
    {
      output << v.S[i] << ", ";
    }
  output << v.S[i] << ")";

  return output;
}
