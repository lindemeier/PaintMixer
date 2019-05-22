#include <PaintMixer/PaintMixer.hxx>
#include <PaintMixer/Palette.hxx>
#include <PaintMixer/Serialization.hxx>

#include <cxxopts.hpp>

#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
  cxxopts::Options options(argv[0], " - Paint mixer command line options");
  options.positional_help("[optional args]").show_positional_help();

  options.add_options()
    // clang-format off
      ("b,basepigments", "path to the the base pigment set file (json)", cxxopts::value<std::string>())     
      ("i,image", "input picture", cxxopts::value<std::string>())
      ("o,output", "Output file to store the extracted palette", cxxopts::value<std::string>()
          ->default_value("extractedPalette.json")->implicit_value("extractedPalette.json"))      
      ("help", "Print help")
      ("n", "desired number of pigments in the extracted palette", cxxopts::value<int>())
      ;
  // clang-format on

  const auto result = options.parse(argc, argv);

  if (result.count("help"))
    {
      std::cout << options.help({"", "Group"}) << std::endl;
      exit(EXIT_SUCCESS);
    }

  if (result.count("b") == 0UL)
    {
      std::cerr << "no base pigments file given" << std::endl;
      std::cout << options.help({"", "Group"}) << std::endl;
      exit(EXIT_FAILURE);
    }

  if (result.count("i") == 0UL)
    {
      std::cerr << "no input picture given" << std::endl;
      std::cout << options.help({"", "Group"}) << std::endl;
      exit(EXIT_FAILURE);
    }
  if (result.count("n") == 0UL)
    {
      std::cerr << "Please define the number of the extracted colors"
                << std::endl;
      std::cout << options.help({"", "Group"}) << std::endl;
      exit(EXIT_FAILURE);
    }

  PaintMixer::Palette basePigments;
  {
    const auto    basePigmentsFile = result["basepigments"].as<std::string>();
    std::ifstream istream(basePigmentsFile);
    try
      {
        PaintMixer::LoadPalette(istream, basePigments);
      }
    catch (cereal::Exception& e)
      {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
      }
  }

  cv::Mat_<PaintMixer::vec3f> image;
  {
    const auto imageFile = result["image"].as<std::string>();
    try
      {
        PaintMixer::LoadImage(imageFile, image);
      }
    catch (std::exception& e)
      {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
      }
  }

  const auto mixer = PaintMixer::PaintMixer(basePigments);

  try
    {
      const auto palette =
        mixer.mixFromInputPicture(image, result["n"].as<int32_t>());

      const auto    outputFile = result["output"].as<std::string>();
      std::ofstream ostream(outputFile);
      PaintMixer::SavePalette(ostream, palette);

      PaintMixer::SaveImage(outputFile + ".basePigments.jpg",
                            VisualizePalette(basePigments, 1.0));

      PaintMixer::SaveImage(outputFile + ".jpg",
                            VisualizePalette(palette, 1.0));
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }

  exit(EXIT_SUCCESS);
}
