#include <random>
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>
#include <omp.h>
#include <vector>

#include "raytracer.hpp"

using namespace LiteMath;
using namespace LiteImage;

constexpr float EPS = 0.001f;
constexpr int max_steps = 16;
std::vector<std::mt19937> generators(omp_get_max_threads());

static float2
mul2x2x2(float2 m[2], const float2 &v)
{
  return m[0]*v[0]+m[1]*v[1];
}

float2 bezier_project(
  const float4 &P1,
  const float4 &P2,
  const float4 &point) {

  return float2 { dot(point, P1), dot(point, P2) };
}

std::optional<float2> trace_surface_newton(
    const float3 &pos,
    const float3 &ray,
    const RBezierGrid &surf) {
  float3 ortho_dir1 = (ray.x || ray.z) ? float3{ 0, 1, 0 } : float3{ 1, 0, 0 };
  float3 ortho_dir2 = normalize(cross(ortho_dir1, ray));
  ortho_dir1 = normalize(cross(ray, ortho_dir2));
  assert(dot(ortho_dir1, ortho_dir2) < 1e-2f);
  assert(dot(ortho_dir1, ray) < 1e-2f);
  assert(dot(ortho_dir2, ray) < 1e-2f);

  float4 P1 = to_float4(ortho_dir1, -dot(ortho_dir1, pos));
  float4 P2 = to_float4(ortho_dir2, -dot(ortho_dir2, pos));
  assert(dot(P1, to_float4(pos, 1.0f)) < 1e-2);
  assert(dot(P2, to_float4(pos, 1.0f)) < 1e-2);
  
  auto distr = std::uniform_real_distribution<float>(0.0f, 1.0f);
  auto &generator = generators[omp_get_thread_num()];
  float2 uv = float2{ distr(generator), distr(generator) };
  float4 surf_point = surf.get_point(uv.x, uv.y);
  surf_point /= surf_point.w;
  float2 D = bezier_project(P1, P2, surf_point);
  
  int steps_left = max_steps-1;
  while(length(D) > EPS && steps_left--) {
    float2 J[2] = 
    { 
      bezier_project(P1, P2, surf.uder(uv.x, uv.y)), //col1
      bezier_project(P1, P2, surf.vder(uv.x, uv.y)) //col2
    };

    float det = J[0][0]*J[1][1] - J[0][1] * J[1][0];

    float2 J_inversed[2] = 
    {
      { J[1][1]/det, -J[0][1]/det },
      { -J[1][0]/det, J[0][0]/det }
    };

    uv = uv - mul2x2x2(J_inversed, D);
    uv.x = clamp(uv.x, 0.0f, 1.0f);
    uv.y = clamp(uv.y, 0.0f, 1.0f);
    assert(0 <= uv.x && uv.x <= 1);
    assert(0 <= uv.y && uv.y <= 1);

    surf_point = surf.get_point(uv.x, uv.y);
    surf_point /= surf_point.w;
    float2 new_D = bezier_project(P1, P2, surf_point);
    
    if (length(new_D) > length(D))
      return {};
    
    D = new_D;
  }

  if (length(D) > EPS)
    return {};
  
  return uv;
}

void draw_points(
    const RBezierGrid &surface,
    const Camera &camera,
    Image2D<uint32_t> &image,
    std::function<ShadeFuncType> shade_function) {
  image.clear(LiteMath::uchar4{ 153, 153, 153, 255 }.u32);
  if (surface.grid.get_n() == 0)
    return;
  float3 forward = normalize(camera.target - camera.position);
  float4x4 view = lookAt(camera.position, camera.target, camera.up);
  float4x4 proj = perspectiveMatrix(camera.fov*180*M_1_PI, camera.aspect, 0.01f, 150.0f);
  int count = 250;
  for (int ui = 0; ui < count; ++ui)
  for (int vi = 0; vi < count; ++vi)
  {
    float u = ui * 1.0f/count;
    float v = vi * 1.0f/count;
    
    float4 point = surface.get_point(u, v);
    point = proj * view * point;
    float w = point.w;
    if (!(-w <= point.x && point.x <= w) ||
        !(-w <= point.y && point.y <= w) ||
        !(-w <= point.z && point.z <= w)) {
      continue;
    }
    point /= point.w;
    float3 normal = surface.normal(u, v);

    uint32_t x = clamp(static_cast<uint32_t>((point.x+1.0f)/2 * image.width()), 0u, image.width()-1);
    uint32_t y = static_cast<uint32_t>((point.y+1.0f)/2 * image.height());
    y = clamp(image.height() - y, 0u, image.height()-1);

    float4 color = shade_function(camera.position, to_float3(point), normal, float2{u, v});
    image[uint2{x, y}] = uchar4{ 
        static_cast<u_char>(color.x*255.0f),
        static_cast<u_char>(color.y*255.0f),
        static_cast<u_char>(color.z*255.0f),
        static_cast<u_char>(1*255.0f) }.u32;
  }
}

void draw_newton(
    const RBezierGrid &surface,
    const Camera &camera,
    Image2D<uint32_t> &image,
    std::function<ShadeFuncType> shade_function) {
  image.clear(LiteMath::uchar4{ 153, 153, 153, 255 }.u32);
  if (surface.grid.get_n() == 0)
    return;
  float4x4 mat  = perspectiveMatrix(camera.fov*180*M_1_PI, camera.aspect, 0.001f, 100.0f)
                * lookAt(camera.position, camera.target, camera.up);
  float4x4 inversed_mat = inverse4x4(mat);
  #pragma omp parallel for 
  for (uint32_t y = 0; y < image.height(); ++y)
  for (uint32_t x = 0; x < image.width();  ++x)
  {
    float2 ndc_point  = float2{ x+0.5f, y+0.5f } 
                      / float2{ image.width()*1.0f, image.height()*1.0f }
                      * 2.0f
                      - 1.0f;
    float4 ndc_point4 = { ndc_point.x, ndc_point.y, 0.0f, 1.0f };
    float4 point = inversed_mat * ndc_point4;
    point /= point.w;

    float3 ray = normalize(to_float3(point)-camera.position);
    float3 pos = camera.position;
    if (!surface.bbox.intersects(pos, ray))
      continue;
    auto intersect_point = trace_surface_newton(pos, ray, surface);
    if (intersect_point.has_value()) {
      float2 uv = intersect_point.value();
      float4 point = surface.get_point(uv.x, uv.y);
      point /= point.w; 
      float3 normal = surface.normal(uv.x, uv.y);
      float4 color = shade_function(camera.position, to_float3(point), normal, uv);
      image[uint2{ x, image.height()-1-y }] = uchar4{ 
        static_cast<u_char>(color[0]*255.0f),
        static_cast<u_char>(color[1]*255.0f),
        static_cast<u_char>(color[2]*255.0f),
        static_cast<u_char>(color[3]*255.0f) }.u32;
    }
  }
}