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

  int width = 500;
  int height = 500;

  Image2D<uint32_t> img(width, height);

  NURBSCurve2D curve = load_nurbs_curve(proj_path / "resources" / "heart.nurbsc");
  auto rbcurves = decompose_nurbs_curve(curve);

  std::vector<RBCurve2DMonotonic> mrbcurves;
  for (auto &curve: rbcurves) {
    mrbcurves.emplace_back(curve, curve.tmin, curve.tmax);
  }

  for (int y = 0; y < height; ++y)
  {
    float u = y * 1.0f / height;
    std::vector<float> vs;
    for (auto &c: mrbcurves) {
      auto cur_vs = c.intersections(u);
      std::copy(cur_vs.begin(), cur_vs.end(), std::back_inserter(vs));
    }
    vs.push_back(0.0f);
    vs.push_back(1.0f);
    std::sort(vs.begin(), vs.end());
    for (int i = 1; i < vs.size(); ++i) {
      if (std::abs(vs[i]-vs[i-1]) < 0.0000001f) {
        vs[i] = vs[i-1];
      }
    }
    vs.resize(std::unique(vs.begin(), vs.end())-vs.begin());
    vs.back() = 1.0f;
    for (int span = 0; span < vs.size()-1; ++span) {
      uchar4 color = (span % 2 == 0) ? uchar4{0, 255, 0, 255} : uchar4{255, 0, 0, 255};
      int x_min = static_cast<int>(vs[span]*width);
      int x_max = static_cast<int>(vs[span+1]*width);
      for (int x = x_min; x < x_max; ++x) {
        img[int2{x, y}] = color.u32;
      }
    }
  }

  SaveBMP("curve_res.bmp", img.data(), width, height);

  return 0;
}