#include <iostream>
#include <filesystem>
#include <LiteMath.h>
#include <Image2d.h>

#include "curve.hpp"
#include "embree_adaptors.hpp"

using namespace LiteMath;
using namespace LiteImage;

int main(int argc, const char **argv) {
  std::filesystem::path exec_path = std::filesystem::canonical(argv[0]);
  std::filesystem::current_path(exec_path.parent_path().parent_path());
  auto proj_path = std::filesystem::current_path(); 

  NURBSCurve2D circle = load_nurbs_curve(proj_path / "resources" / "circle.nurbsc");
  NURBSCurve2D heart = load_nurbs_curve(proj_path / "resources" / "heart.nurbsc");
  auto rbeziers = circle.decompose();
  auto heart_rbeziers = heart.decompose();
  //rbeziers = heart_rbeziers;
  std::copy(heart_rbeziers.begin(), heart_rbeziers.end(), std::back_inserter(rbeziers));

  auto [boxes, leaves] = get_kdtree_leaves(rbeziers);
  std::cout << boxes.size() << std::endl;

  int w = 1000, h = 1000;
  Image2D<uint32_t> img(w, h), img2(w, h);


  embree::EmbreeTrimKdTree tree;
  for (int i = 0; i < boxes.size(); ++i) {
    tree.add_box(rbeziers.data(), boxes[i], leaves[i]);
  }
  tree.commit_scene();

  for (int y = 0; y < w; ++y)
  for (int x = 0; x < h; ++x)
  {
    float u = (y+0.5f)/h;
    float v = (x+0.5f)/w;

    bool inside = tree.query(u, v);
    if (inside) {
      img[int2{x, y}] = 0xff0000ff;
    } else {
      img[int2{x, y}] = 0xff00ff00;
    }
  }

  auto save_path = proj_path / "result.bmp";
  SaveBMP(save_path.c_str(), img.data(), w, h);
  return 0;
}
