#include <cstddef>

#include "curve.hpp"
#include "Image2d.h"
#include "LiteMath.h"

using namespace LiteMath;
using namespace LiteImage;

int main(int argc, const char **argv) {
  std::filesystem::path exec_path = std::filesystem::canonical(argv[0]);
  std::filesystem::current_path(exec_path.parent_path().parent_path());
  auto proj_path = std::filesystem::current_path();

  int width = 800;
  int height = 600;

  Image2D<uint32_t> img(width, height);

  NURBSCurve2D curve = load_nurbs_curve(proj_path / "resources" / "circle.nurbsc");
  auto rbcurves = decompose_nurbs_curve(curve);

  std::vector<RBCurve2DMonotonic> mrbcurves;
  for (auto &curve: rbcurves) {
    mrbcurves.emplace_back(curve, curve.tmin, curve.tmax);
  }


  return 0;
}