#include <vector>
#include <cassert>

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

float point(const LiteMath::float3 *pw, size_t index, size_t p, float u)
{
  float u_n = 1.0f;
  float _1_u = 1.0f - u;
  int bc = 1.0f;

  float cw = pw[0][index] * _1_u;
  for (int i = 1; i <= p-1; ++i) {
    u_n *= u;
    bc = bc * (p-i+1)/i;
    cw = (cw + u_n * bc * pw[i][index]) * _1_u;
  }
  cw += (u_n * u) * pw[p][index];
  return cw;
}

void __enzyme_autodiff(...);

template<typename RT, typename... Args>
RT __enzyme_autodiff(void*, Args...);

int enzyme_dup;
int enzyme_out;
int enzyme_const;

float der(const LiteMath::float3 *pw, size_t index, size_t p, float u)
{
  float res = __enzyme_autodiff<float>(
    (void*)point,
    enzyme_const, pw, 
    enzyme_const, index,
    enzyme_const, p,
    enzyme_out, u);
  return res;
}

float dfg_fdg(const LiteMath::float3 *pw, size_t index, size_t p, float u) {
  float g = point(pw, 2, p, u);
  float dg = der(pw, 2, p, u);
  float f = point(pw, index, p, u);
  float df = der(pw, index, p, u);
  return df * g - f * dg;
}

float der_dfg_fdg(const LiteMath::float3 *pw, size_t index, size_t p, float u) {
  float res = __enzyme_autodiff<float>(
    (void*)dfg_fdg,
    enzyme_const, pw, 
    enzyme_const, index,
    enzyme_const, p,
    enzyme_out, u);
  return res;
}