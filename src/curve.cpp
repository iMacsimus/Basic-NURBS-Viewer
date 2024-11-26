#include <vector>
#include <cassert>
#include <optional>

#include "curve.hpp"

using namespace LiteMath;

int RBCurve2D::degree() const {
  return pw.size()-1;
}

LiteMath::float3 RBCurve2D::get_point(float u) const {
  int p = degree();

  float u_n = 1.0f;
  float _1_u = 1.0f - u;
  int bc = 1.0f;

  float3 res = pw[0] * _1_u;
  for (int i = 1; i <= p-1; ++i) {
    u_n *= u;
    bc = bc * (p-i+1)/i;
    res = (res + u_n * bc * pw[i]) * _1_u;
  }
  res += (u_n * u) * pw[p];
  return res;
}

LiteMath::float3 RBCurve2D::non_rat_der(float u) const {
  int p = degree();

  if (p == 1) {
    return p * (pw[1]-pw[0]);
  }

  float u_n = 1.0f;
  float _1_u = 1.0f - u;
  int bc = 1.0f;
  
  float3 next = pw[1];
  float3 cur = pw[0];
  float3 res = (next-cur) * _1_u;
  cur = next;

  for (int i = 1; i <= p-2; ++i) {
    u_n *= u;
    bc = bc * (p-i)/i;
    next = pw[i+1];
    res = (res + u_n * bc * (next-cur)) * _1_u;
    cur = next;
  }

  next = pw[p];
  res += (u_n * u) * (next-cur);

  res *= p;
  return res;
}

LiteMath::float3 RBCurve2D::der(float u) const {
  float3 cw = get_point(u);
  float3 nr_der = non_rat_der(u);
  return (nr_der * cw.z - cw * nr_der.z) / (cw.z * cw.z);
}

Matrix2D<float3>
nr_derivative_points(const RBCurve2D &curve)
{
  int p = curve.degree();
  Matrix2D<float3> res(p+1, p+1);
  std::copy(curve.pw.begin(), curve.pw.end(), res.data());

  for (int order = 1; order <= p; ++order) {
    for (int i = 0; i <= p-order; ++i) {
      res[{order, i}] = (p-order+1) * (res[{order-1, i+1}] - res[{order-1, i}]);
    }
  }

  return res;
}

float3 get_nr_der_point(const Matrix2D<float3> &ders_points, int order, float u)
{
  int p = ders_points.get_n()-1;
  if (order == p) {
    return ders_points[{order, 0}];
  }
  if (order > p) {
    return float3{ 0.0f };
  }

  const float3 *pw = &ders_points[{order, 0}];
  p -= order;

  float u_n = 1.0f;
  float _1_u = 1.0f - u;
  int bc = 1.0f;

  float3 res = pw[0] * _1_u;
  for (int i = 1; i <= p-1; ++i) {
    u_n *= u;
    bc = bc * (p-i+1)/i;
    res = (res + u_n * bc * pw[i]) * _1_u;
  }
  res += (u_n * u) * pw[p];
  return res;
}

float3 fder_g(const Matrix2D<float3> &der_points, int order, float u) {
  int bc = 1;
  float3 res = {};
  for (int i = 0; i <= order; ++i) {
    float3 f_der_i = get_nr_der_point(der_points, i+1, u);
    float g_n_i = get_nr_der_point(der_points, order-i, u).z;
    res += bc * g_n_i * f_der_i;
    bc = bc * (order-i) / (i+1);
  }
  return res;
}

float3 f_gder(const Matrix2D<float3> &der_points, int order, float u) {
  int bc = 1;
  float3 res = {};
  for (int i = 0; i <= order; ++i) {
    float g_der_i = get_nr_der_point(der_points, i+1, u).z;
    float3 f_n_i = get_nr_der_point(der_points, order-i, u);
    res += bc * g_der_i * f_n_i;
    bc = bc * (order-i) / (i+1);
  }
  return res;
}

float3 fder_g__f_gder(const Matrix2D<float3> &der_points, int order, float u) {
  float3 a = fder_g(der_points, order, u);
  float3 b = f_gder(der_points, order, u);
  return a-b;
}

template<typename F>
std::optional<float>
bisection(F f, float min_u, float max_u) {
  float l = min_u;
  float r = max_u;
  {
    float tmp1 = f(min_u);
    float tmp2 = f(max_u);
    if (tmp1 >= 0 && tmp2 >= 0) {
      return {};
    }
    if (tmp1 <= 0 && tmp2 <= 0) {
      return {};
    }
    if (tmp1 > tmp2) {
      std::swap(l, r);
    }
  }

  constexpr float BISECTION_EPS = 0.001f;
  while(std::abs(l-r) < BISECTION_EPS) {
    float m = (l+r) / 2.0f;
    float m_val = f(m);
    if (m_val < 0.0f) {
      l = m;
    } else {
      r = m;
    }
  }

  return (l+r) / 2.0f;
}

std::vector<float>
get_monotonic_parts(const Matrix2D<float3> &der_points, int index, int order) {
  int p = der_points.get_n() - 1;
  if (order == 2*p-1) { 
    return { 0.0f, 1.0f };
  }

  std::vector<float> res = { 0.0f };

  auto f = [&](float u) {
    return fder_g__f_gder(der_points, order, u)[index];
  };

  auto knots = get_monotonic_parts(der_points, index, order+1);
  for (int span = 0; span < knots.size()-1; ++span) {
    auto potential_root = bisection(f, knots[span], knots[span+1]);
    if (potential_root.has_value()) {
      res.push_back(potential_root.value());
    }
  }
  
  res.push_back(1.0f);
  return res;
}

std::vector<float>
get_bimonotonic_parts(const RBCurve2D &curve) {
  auto der_matrix = nr_derivative_points(curve);
  auto knots1 = get_monotonic_parts(der_matrix, 0, 0);
  auto knots2 = get_monotonic_parts(der_matrix, 1, 0);
  std::vector<float> knots;
  std::merge(knots1.begin(), knots1.end(), knots2.begin(), knots2.end(), std::back_inserter(knots));
  knots.resize(std::unique(knots.begin(), knots.end())-knots.begin());
  return knots;
}

