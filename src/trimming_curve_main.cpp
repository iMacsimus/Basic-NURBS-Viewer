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

  int points_count = 1000;
  bool even = true;
  for (auto &c: mrbcurves) {
    for (int span = 0; span < c.monotonic_knots.size()-1; ++span) {
      float tmin = lerp(c.tmin, c.tmax, c.monotonic_knots[span]);
      float tmax = lerp(c.tmin, c.tmax, c.monotonic_knots[span+1]);
      for (int point_i = 0; point_i < points_count; ++point_i) {
        float t = lerp(tmin, tmax, point_i * 1.0f/points_count);
        auto p = c.get_point(t);
        p /= p.z;
        int x = static_cast<int>(p.y*width);
        int y = static_cast<int>(p.x*height);
        x = clamp(x, 0, width-1);
        y = clamp(y, 0, width-1);
        for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx)
        {
          int2 xy = { x+dx, y+dy };
          if (any_of(xy < int2{0, 0}) || any_of(xy >= int2{width, height})) {
            continue;
          }
          uchar4 color = even ? uchar4{0, 0, 0, 255} : uchar4{255, 255, 255, 255};
          img[xy] = color.u32;
        }
      }
      even = !even;
    }
  }

  SaveBMP("curve_res.bmp", img.data(), width, height);

  return 0;
}