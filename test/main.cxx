#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include <PaintMixer.hxx>
#include <Palette.hxx>
#include <Serialization.hxx>

#include <iostream>
#include <sstream>

PaintMixer::Palette CreateCurtisPigments()
{ /*
   * Cassidy J. Curtis, Sean E. Anderson, Joshua E. Seims, Kurt W.
   * Fleischer, and David H. Salesin. 1997. Computer-generated
   * watercolor. In Proceedings of the 24th annual conference on
   * Computer graphics and interactive techniques (SIGGRAPH '97). ACM
   * Press/Addison-Wesley Publishing Co., New York, NY, USA, 421-430.
   * DOI=http://dx.doi.org/10.1145/258734.258896
   * */
  //    http://grail.cs.washington.edu/projects/watercolor/paper_small.pdf

  PaintMixer::Palette palette;
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.0, 0.5, 1.0}, {0.2, 0.3, 0.4}});

  //    “Quinacridone Rose”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.22, 1.47, 0.57}, {0.05, 0.003, 0.03}});

  //    “Indian Red”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.46, 1.07, 1.50}, {1.28, 0.38, 0.21}});

  //    “Cadmium Yellow”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.10, 0.36, 3.45}, {0.97, 0.65, 0.007}});

  //    “Hookers Green”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{1.62, 0.61, 1.64}, {0.01, 0.012, 0.003}});

  //    “Cerulean Blue”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{1.52, 0.32, 0.25}, {0.06, 0.26, 0.40}});

  //    “Burnt Umber”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.74, 1.54, 2.10}, {0.09, 0.09, 0.004}});

  //    “Cadmium Red”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.14, 1.08, 1.68}, {0.77, 0.015, 0.018}});

  //    “Brilliant Orange”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.13, 0.81, 3.45}, {0.005, 0.009, 0.007}});

  //    “Hansa Yellow”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.06, 0.21, 1.78}, {0.50, 0.88, 0.009}});

  //    “Phthalo Green”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{1.55, 0.47, 0.63}, {0.01, 0.05, 0.035}});

  //    “French Ultramarine”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.86, 0.86, 0.06}, {0.005, 0.005, 0.09}});

  //    “Interference Lilac”
  palette.emplace_back(
    PaintMixer::PaintCoeff{{0.08, 0.11, 0.07}, {1.25, 0.42, 1.43}});

  return palette;
};

TEST_CASE("Palette serialization")
{
  PaintMixer::Palette palette = CreateCurtisPigments();

  std::stringstream out;
  PaintMixer::SavePalette(out, palette);

  std::stringstream   in(out.str());
  PaintMixer::Palette loadedPalette;
  PaintMixer::LoadPalette(in, loadedPalette);

  REQUIRE(palette.size() == loadedPalette.size());
  for (auto l = 0U; l < palette.size(); l++)
    {
      for (auto i = 0U; i < PaintMixer::CoeffSamplesCount; i++)
        {
          REQUIRE(palette[l].K[i] == Approx(loadedPalette[l].K[i]));
          REQUIRE(palette[l].S[i] == Approx(loadedPalette[l].S[i]));
        }
    }
}

TEST_CASE("Paint mixing")
{
  const auto testPaint = [](const PaintMixer::PaintCoeff& expectedPaint,
                            const PaintMixer::PaintCoeff& mixedPaint) {
    for (auto i = 0U; i < PaintMixer::CoeffSamplesCount; i++)
      {
        REQUIRE(expectedPaint.K[i] == Approx(mixedPaint.K[i]));
        REQUIRE(expectedPaint.S[i] == Approx(mixedPaint.S[i]));
      }
  };

  const auto basePalette = CreateCurtisPigments();
  auto       mixer       = PaintMixer::PaintMixer(basePalette);

  // test getter
  {
    const auto palette = mixer.getUnderlyingPalette();
    REQUIRE(palette.size() == basePalette.size());
    for (auto l = 0U; l < palette.size(); l++)
      {
        testPaint(palette[l], basePalette[l]);
      }
  }

  // mix a base paint
  {
    const auto baseCount     = mixer.getUnderlyingPalette().size();
    const auto expectedPaint = mixer.getUnderlyingPalette().front();

    {
      std::vector<float64_t> weights(baseCount, 0.0);
      weights.front() = 1.0;

      const auto mixedPaint = mixer.mixSinglePaint(weights);
      testPaint(expectedPaint, mixedPaint);
    }
    {
      std::vector<float64_t> weights(baseCount, 0.0);
      weights.front() = 0.5;

      const auto mixedPaint = mixer.mixSinglePaint(weights);
      testPaint(expectedPaint, mixedPaint);
    }
  }

  // mixing get the mixing recipe for a previously mixed paint
  {
    const auto             baseCount = mixer.getUnderlyingPalette().size();
    std::vector<float64_t> expectedWeights(baseCount, 0.0);
    expectedWeights.front() = 0.2;
    expectedWeights[4U]     = 1.0 - expectedWeights.front();

    const auto mixedPaint = mixer.mixSinglePaint(expectedWeights);

    const auto recipe = mixer.getWeightsForMixingTargetPaint(mixedPaint);
  }
}