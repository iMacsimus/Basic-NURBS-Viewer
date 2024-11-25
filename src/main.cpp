#include <iostream>
#include <filesystem>
#include <cmath>
#include <cinttypes>
#include <string_view>
#include <cstdio>
#include <random>

#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"

#include <SDL2/SDL.h>
#include <Image2d.h>
#include <LiteMath.h>

#include "Surface.hpp"
#include "raytracer.hpp"
#include "utils.hpp"
#include "embree_adaptors.hpp"
#include "gui.hpp"
#include "curve.hpp"

using namespace LiteMath;
using namespace LiteImage;

int main(int, char** argv)
{
  std::filesystem::path exec_path = std::filesystem::canonical(argv[0]);
  std::filesystem::current_path(exec_path.parent_path().parent_path());
  auto proj_path = std::filesystem::current_path();

  RBCurve2D rbcurve;
  rbcurve.pw = {
    { -1, -1, 1},
    { -1, 1, 1},
    { 1, 1, 1},
    { 1, -1, 1}
  };
  srand(time(NULL));
  float u = rand() * 1.0f / RAND_MAX;
  std::cout << der(rbcurve.pw.data(), 0, 3, u) << " ";
  std::cout << der(rbcurve.pw.data(), 1, 3, u) << " ";
  std::cout << der(rbcurve.pw.data(), 2, 3, u) << std::endl;
  float3 der = rbcurve.non_rat_der(u);
  std::cout << der.x << " " << der.y << " " << der.z << std::endl;
}
