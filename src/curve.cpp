#include <vector>
#include <cassert>
#include <optional>
#include <fstream>
#include <iostream>

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

std::vector<RBCurve2DMapped>
decompose_nurbs_curve(const NURBSCurve2D &nurbs) {
  int n = nurbs.pw.size()-1;
  int m = nurbs.knots.size()-1;
  int p = m-n-1;

  const float3 *Pw = nurbs.pw.data();
  const float *U = nurbs.knots.data();

  int a = p;
  int b = p+1;
  int nb = 0;
  std::vector<RBCurve2DMapped> Qw(1, RBCurve2D{std::vector<float3>(p+1)});
  std::vector<float> alphas(p+1);
  for (int i=0; i <= p; ++i) 
    Qw[nb].pw[i] = Pw[i];
  while (b < m) {
    int i = b;
    while (b < m && U[b+1] == U[b]) 
      b++;
    int mult = b-i+1;
    if (mult < p) {
      float numer = U[b]-U[a]; /* Numerator of alpha */
      /* Compute and store alphas */
      for (int j = p; j > mult; j--)
        alphas [j-mult-1] = numer/(U[a+j] - U[a]);
      int r = p-mult; /* Insert knot r times */
      for (int j= 1; j <= r; j++) {
        int save = r-j;
        int s = mult+j; /* This many new points */
        for (int k=p; k>=s; k--) {
          float alpha = alphas[k-s];
          Qw[nb].pw[k] = alpha*Qw[nb].pw[k] + (1.0f - alpha)*Qw[nb].pw[k - 1];
        }
        if (b < m) { /* Control point of */ 
          if (Qw.size() == nb+1)
            Qw.push_back(RBCurve2D{std::vector<float3>(p+1)});
          Qw [nb+1].pw[save] = Qw[nb].pw[p]; /* next segment */
        }
      }
    }
    nb = nb+1; /* Bezier segment completed */
    if (b < m) { 
      if (Qw.size() == nb)
        Qw.push_back(RBCurve2D{std::vector<float3>(p+1)});
      /* Initialize for next segment */
      for (i=p-mult; i<=p; i++) 
        Qw[nb].pw[i] = Pw[b-p+i];
      a = b;
      b = b+1;
    }
  }

  auto knots = nurbs.knots;
  knots.resize(std::unique(knots.begin(),knots.end())-knots.begin());

  for (int i = 0; i < knots.size()-1; ++i) {
    Qw[i].tmin = knots[i];
    Qw[i].tmax = knots[i+1];
  }
  return Qw;
}

NURBSCurve2D load_nurbs_curve(std::filesystem::path path) {
  std::fstream fin;
  fin.exceptions(std::ios::failbit|std::ios::badbit);
  fin.open(path);
  if (!fin.good()) throw std::runtime_error("Failed to open the file.");

  std::string tmp_str;
  char tmp_chr;

  NURBSCurve2D curve;

  int n;
  fin >> tmp_chr >> tmp_chr >> n; // n = ...

  curve.pw = std::vector<float3>(n+1, { 0, 0, 1.0f });
  float min_u = std::numeric_limits<float>::infinity();
  float max_u = -min_u;
  float min_v = min_u;
  float max_v = max_u;

  fin >> tmp_str; // "points:"
  for (int i = 0; i <= n; ++i)
  {
    auto &point = curve.pw[i];
    fin >> tmp_chr >> point.x >> tmp_chr >> point.y >> tmp_chr; // { ..., ..., ... }
    min_u = std::min(min_u, point.x);
    max_u = std::max(max_u, point.x);
    min_v = std::min(min_v, point.y);
    max_v = std::max(max_v, point.y);
  }
  // std::cout << "normalized points:" << std::endl;
  // for (int i = 0; i <= n; ++i)
  // {
  //   auto &point = curve.pw[i];
  //   point.x = (point.x-min_u) / (max_u-min_u);
  //   point.y = (point.y-min_v) / (max_v-min_v);
  //   std::cout << "{" << point.x << ", " << point.y << "} ";
  // }
  // std::cout << std::endl;


  fin >> tmp_str; // "weights:"
  for (int i = 0; i <= n; ++i)
  {
    float w;
    fin >> w;
    curve.pw[i] *= w;
  }

  fin >> tmp_str; // "degree:"
  int deg;
  fin >> deg;

  curve.knots.resize(n+deg+2);
  fin >> tmp_str; // "knots:"
  for (size_t i = 0; i < curve.knots.size(); ++i) {
    fin >> curve.knots[i];
  }

  float u_min = curve.knots.front();
  float u_max = curve.knots.back();
  for (auto &elem: curve.knots)
    elem = (elem-u_min)/(u_max-u_min);

  return curve;
}