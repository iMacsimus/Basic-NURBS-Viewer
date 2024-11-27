#ifndef NURBS_SRC_CURVE
#define NURBS_SRC_CURVE

#include <vector>
#include <filesystem>

#include "utils.hpp"
#include "LiteMath.h"

struct RBCurve2D
{
public:
  int degree() const;
public:
  LiteMath::float3 get_point(float t) const;
  LiteMath::float3 der(float t) const;
  LiteMath::float3 non_rat_der(float t) const;
public:
  std::vector<LiteMath::float3> pw;
};
std::vector<float>
get_bimonotonic_parts(const RBCurve2D &curve);

struct RBCurve2DMapped: RBCurve2D
{
public:
  RBCurve2DMapped() = default;
  RBCurve2DMapped(const RBCurve2D &curve, float tmin = 0.0f, float tmax = 1.0f)
      : RBCurve2D(curve), tmin(tmin), tmax(tmax) {}
public:
  float tmin = 0.0f; 
  float tmax = 1.0f;
public:
  LiteMath::float3 get_point(float t) const {
    t = (t - tmin) / (tmax-tmin);
    return RBCurve2D::get_point(t);
  }
  LiteMath::float3 der(float t) const {
    t = (t - tmin) / (tmax-tmin);
    return RBCurve2D::der(t);
  }
  LiteMath::float3 non_rat_der(float t) const {
    t = (t - tmin) / (tmax-tmin);
    return RBCurve2D::non_rat_der(t);
  }
};

struct RBCurve2DMonotonic: RBCurve2DMapped
{
public:
  RBCurve2DMonotonic() = default;
  RBCurve2DMonotonic(const RBCurve2D &curve, float tmin = 0.0f, float tmax = 1.0f)
      : RBCurve2DMapped(curve, tmin, tmax)
      , monotonic_knots(get_bimonotonic_parts(curve)) {}
public:
  std::vector<float> monotonic_knots;
};

struct NURBSCurve2D
{
public:
  int degree() const { return knots.size()-pw.size()-1; }
  int get_n() const { return pw.size()-1; }
public:
  std::vector<LiteMath::float3> pw;
  std::vector<float> knots;
};
NURBSCurve2D load_nurbs_curve(std::filesystem::path path);
std::vector<RBCurve2DMapped> decompose_nurbs_curve(const NURBSCurve2D &nurbs);

#endif 