#ifndef NURBS_SRC_SURFACE
#define NURBS_SRC_SURFACE

#include <cinttypes>
#include <vector>
#include <cmath>
#include <filesystem>
#include <limits>
#include <iostream>

#include <LiteMath.h>
#include <Image2d.h>

#include <LiteMath.h>
#include <Image2d.h>

template<typename T>
struct Matrix2D 
{
public:
  Matrix2D(int n = 0, int m = 0, T value = T()): n(n), m(m), values(n * m, value) {}
  T& operator[](std::pair<int, int> ids) { return values[ids.first*m+ids.second]; }
  const T& operator[](std::pair<int, int> ids) const { return values[ids.first*m+ids.second]; }
  const T* data() const { return values.data(); }
  int get_n() const { return n; }
  int get_m() const { return m; }
private:
  int n; int m;
  std::vector<T> values;
};

struct BoundingBox3d
{
public:
  BoundingBox3d() = default;
  BoundingBox3d(const LiteMath::float4 *points, int count) {
    for (int i = 0; i < count; ++i) {
      LiteMath::float3 point = LiteMath::to_float3(points[i]/points[i].w);
      mn = LiteMath::min(mn, point);
      mx = LiteMath::max(mx, point);
    }
  }
public:
  LiteMath::float3 mn = LiteMath::float3{std::numeric_limits<float>::infinity()};
  LiteMath::float3 mx = -LiteMath::float3{std::numeric_limits<float>::infinity()};
  LiteMath::float3 center() const { return (mn+mx)/2; }
  LiteMath::float3 shape() const { return (mn-mx); }
public:
  bool intersects(
      const LiteMath::float3 &pos,
      const LiteMath::float3 &ray) const {
    LiteMath::float3 lo = (mn-pos)/ray;
    LiteMath::float3 ho = (mx-pos)/ray;
    LiteMath::float3 tmin3 = LiteMath::min(lo, ho);
    LiteMath::float3 tmax3 = LiteMath::max(lo, ho);
    float tmin = LiteMath::hmax(tmin3);
    float tmax = LiteMath::hmin(tmax3);
    if (tmin > tmax)
      return false;
    return tmin >= 0.0f;
  }
};

class Surface {
public:
  Surface() = default;
  Surface(
      const Matrix2D<LiteMath::float4> &points,
      const Matrix2D<float> &weights,
      uint32_t deg_u, uint32_t deg_v,
      const std::vector<float> &u_knots,
      const std::vector<float> &v_knots)
      : points(points), weights(weights),
          deg_u(deg_u), deg_v(deg_v),
          u_knots(u_knots), v_knots(v_knots),
          bbox(points.data(), points.get_n()*points.get_m()){}
public:
  Matrix2D<LiteMath::float4> points;
  Matrix2D<float> weights;
  uint32_t deg_u;
  uint32_t deg_v;
  std::vector<float> u_knots;
  std::vector<float> v_knots;
  BoundingBox3d bbox;
};

struct SurfaceView
{
public:
  SurfaceView() = delete;
  SurfaceView(const Surface &surf): p_surf(&surf) {
    left.resize(std::max(p(), q())+1);
    right.resize(std::max(p(), q())+1);
    Nu.resize(p()+1);
    Nv.resize(q()+1);
    temp.resize(q()+1);
  }
public:
  const Surface *p_surf;
  mutable std::vector<float> left, right; //temporary buffers
  mutable std::vector<float> Nu, Nv; //temporary buffers
  mutable std::vector<LiteMath::float4> temp; //temporary buffer
public:
  int n() const { return p_surf->points.get_n()-1; }
  int m() const { return p_surf->points.get_m()-1; }
  int p() const { return p_surf->deg_u; }
  int q() const { return p_surf->deg_v; }
public:
  int find_span(int n, int p, float u, const float *U) const;
  void basis_funs(int i, float u, int p, const float *U, float *N) const;
public:
  LiteMath::float4 get_point(float u, float v) const;
  LiteMath::float4 uderivative(float u, float v) const;
  LiteMath::float4 vderivative(float u, float v) const;
  LiteMath::float3 get_normal(float u, float v) const;
  bool u_closed() const;
  bool v_closed() const;
};

Surface load_surface(const std::filesystem::path &path);

#endif