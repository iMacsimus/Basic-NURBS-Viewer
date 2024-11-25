#ifndef NURBS_SRC_CURVE
#define NURBS_SRC_CURVE

#include <vector>

#include "utils.hpp"
#include "LiteMath.h"

struct RBCurve2D
{
public:
  int degree() const;
public:
  LiteMath::float3 get_point(float u) const;
  LiteMath::float3 non_rat_der(float u) const;
public:
  std::vector<LiteMath::float3> pw;
};

float point(const LiteMath::float3 *pw, size_t index, size_t p, float u);
float der(const LiteMath::float3 *pw, size_t index, size_t p, float u);

float dfg_fdg(const LiteMath::float3 *pw, size_t index, size_t p, float u);
float der_dfg_fdg(const LiteMath::float3 *pw, size_t index, size_t p, float u);

#endif 