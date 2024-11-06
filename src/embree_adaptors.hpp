// Copyright 2009-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
#ifndef NURBS_SRC_EMBREE_ADAPTORS
#define NURBS_SRC_EMBREE_ADAPTORS

#include <list>
#include <embree4/rtcore.h>
#include <LiteMath.h>
#include <Image2d.h>

#include "Surface.hpp"
#include "raytracer.hpp"

namespace embree
{
  struct RBGridView
  {
    const RBezierGrid *p_grid;
    const std::vector<BoundingBox3d> *p_boxes;
    const std::vector<LiteMath::float2> *p_uvs;
  };
  void rbgrid_bounds_function(const RTCBoundsFunctionArguments *args);
  void rbgrid_intersect_function(const RTCIntersectFunctionNArguments *args);
  void rbgrid_occluded_function(const RTCOccludedFunctionNArguments *args);

  void errorFunction(void* userPtr, enum RTCError error, const char* str);
  struct EmbreeScene
  {
  public:
    EmbreeScene() {
      device = rtcNewDevice(nullptr);
#ifndef NDEBUG
      rtcSetDeviceErrorFunction(device, errorFunction, nullptr);
#endif
      scn = rtcNewScene(device);
    }
  public:
    void attach_surface(
        const RBezierGrid &rbezier, 
        const std::vector<BoundingBox3d> &boxes,
        const std::vector<LiteMath::float2> &uvs) {
      views.push_back(RBGridView{ &rbezier, &boxes, &uvs });
      auto &view = views.back();
      RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
      rtcSetGeometryUserPrimitiveCount(geom, boxes.size());
      rtcSetGeometryUserData(geom, &view);
      rtcSetGeometryBoundsFunction(geom, rbgrid_bounds_function, nullptr);
      rtcSetGeometryIntersectFunction(geom, rbgrid_intersect_function);
      //rtcSetGeometryOccludedFunction(geom, rbgrid_occluded_function);
      rtcCommitGeometry(geom);
      rtcAttachGeometry(scn, geom);
      rtcReleaseGeometry(geom);
    }
  public:
    void clear_scene() {
      rtcReleaseScene(scn);
      views.resize(0);
      scn = rtcNewScene(device);
    }
  public:
    void commit_scene() {
      rtcCommitScene(scn);
    }
  public:
    ~EmbreeScene() {
      rtcReleaseScene(scn);
      rtcReleaseDevice(device);
    }
  public:
    void draw(const Camera &camera, FrameBuffer &fb, std::function<ShadeFuncType> shade_func = shade_uv) const;
  private:
    std::list<RBGridView> views;
    RTCDevice device;
    RTCScene scn;
  };
} // namespace embree

#endif

