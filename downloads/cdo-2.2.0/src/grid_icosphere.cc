// modified code from:
// https://schneide.wordpress.com/2016/07/15/generating-an-icosphere-in-c

#include "vector3d.h"

#include <utility>
#include <limits>
#include <cstdio>
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>

using Index = size_t;
using Vertex = Vector3d;
using Triangle = std::array<Index, 3>;
using TriangleList = std::vector<Triangle>;
using VertexList = std::vector<Vertex>;

using Lookup = std::map<std::pair<Index, Index>, Index>;

// icosahedron from ICON
namespace icosahedron
{

// Northern hemisphere are the first 6 elements of vertices[0:5]
// Southern hemisphere are the other 6 elements of vertices[6:11]
// 12 vertices
static VertexList vertices(12);

void
init(void)
{
  constexpr auto pi_5 = M_PI * 0.2;
  // first define the vertices of the icosahedron
  const auto z_w = 2.0 * std::acos(1.0 / (2.0 * std::sin(pi_5)));

  // set poles first - it is simple
  vertices[0] = Vertex{ 0.0, 0.0, 1.0 };
  vertices[11] = Vertex{ 0.0, 0.0, -1.0 };
  // now set the vertices on the two latitude rings
  int i_mdist[10];
  for (int j = 1; j < 11; ++j)
    {
      if (j % 2 == 0)
        i_mdist[j / 2 + 4] = -1 + (j - 1) - 10 * ((j - 1) / 7);
      else
        i_mdist[(j + 1) / 2 - 1] = -1 + (j - 1) - 10 * ((j - 1) / 7);
    }

  for (int j = 1; j < 11; ++j)
    {
      // toggle the hemisphere
      const auto i_msgn = (j >= 6) ? -1.0 : 1.0;
      // compute the meridian angle for the base vertex.
      const auto z_rlon = (1.0 + i_mdist[j - 1]) * pi_5;
      // now initialize the coordinates
      vertices[j] = Vertex{ std::sin(z_w) * std::cos(z_rlon), std::sin(z_w) * std::sin(z_rlon), std::cos(z_w) * i_msgn };
    }
}

// 20 triangles
static const TriangleList triangles
    = { { { 0, 1, 2 } },  { { 0, 2, 3 } },  { { 0, 3, 4 } },  { { 0, 4, 5 } },   { { 0, 5, 1 } },
        { { 6, 2, 1 } },  { { 7, 3, 2 } },  { { 8, 4, 3 } },  { { 9, 5, 4 } },   { { 10, 1, 5 } },
        { { 2, 6, 7 } },  { { 3, 7, 8 } },  { { 4, 8, 9 } },  { { 5, 9, 10 } },  { { 1, 10, 6 } },
        { { 11, 7, 6 } }, { { 11, 8, 7 } }, { { 11, 9, 8 } }, { { 11, 10, 9 } }, { { 11, 6, 10 } } };

}  // namespace icosahedron

static Index
vertexForEdge(Lookup &lookup, VertexList &vertices, Index first, Index second)
{
  Lookup::key_type key(first, second);
  if (key.first > key.second) std::swap(key.first, key.second);

  const auto inserted = lookup.insert({ key, vertices.size() });
  if (inserted.second) vertices.push_back((vertices[first] + vertices[second]).normalised());

  return inserted.first->second;
}

static TriangleList
subdivide(VertexList &vertices, const TriangleList &triangles)
{
  Lookup lookup;
  Triangle mid;

  const auto n = triangles.size();
  TriangleList result(4 * n);
  for (size_t i = 0; i < n; ++i)
    {
      const auto &each = triangles[i];
      for (int edge = 0; edge < 3; edge++) { mid[edge] = vertexForEdge(lookup, vertices, each[edge], each[(edge + 1) % 3]); }

      result[i * 4 + 0] = { each[0], mid[0], mid[2] };
      result[i * 4 + 1] = { each[1], mid[1], mid[0] };
      result[i * 4 + 2] = { each[2], mid[2], mid[1] };
      result[i * 4 + 3] = mid;
    }

  return result;
}

using IndexedMesh = std::pair<VertexList, TriangleList>;

static IndexedMesh
makeIcosphere(int subdivisions)
{
  icosahedron::init();
  auto vertices = icosahedron::vertices;
  auto triangles = icosahedron::triangles;

  for (int i = 0; i < subdivisions; ++i) triangles = subdivide(vertices, triangles);

  return { vertices, triangles };
}

static inline void
d_normalize(Vertex &v)
{
  const auto dnorm = std::sqrt(v * v);
  v /= dnorm;
}

static Vertex
circumCenterMean(const Vertex &v0, const Vertex &v1, const Vertex &v2)
{
  /*
    v0, v1, v2: the coordinates of the three triangle vertices (_dmo,nit vectors) in
    counter clockwise order center: the coordinates of the circumcenter unless co-linear
  */
  // cu0, cu1, cu2: vector product of center:  e1 x e2
  // e1, e2: edges of the underlying planar triangle: v1-v0 ands v2-v0, respectively

  auto e1 = v1 - v0;
  auto e2 = v2 - v0;
  auto cu0 = e1 % e2;
  if ((cu0 * v0) < 0.0) cu0 = -cu0;
  d_normalize(cu0);

  e1 = v2 - v1;
  e2 = v0 - v1;
  auto cu1 = e1 % e2;
  if ((cu1 * v1) < 0.0) cu1 = -cu1;
  d_normalize(cu1);

  e1 = v0 - v2;
  e2 = v1 - v2;
  auto cu2 = e1 % e2;
  if ((cu2 * v2) < 0.0) cu2 = -cu2;
  d_normalize(cu2);

  auto center = cu0 + cu1 + cu2;
  d_normalize(center);

  return center;
}

size_t
genIcosphereCoords(int subdivisions, bool withBounds, std::vector<double> &xvals, std::vector<double> &yvals,
                   std::vector<double> &xbounds, std::vector<double> &ybounds)
{
  const auto mesh = makeIcosphere(subdivisions);
  const auto &vertices = mesh.first;
  const auto &triangles = mesh.second;

  const auto ncells = triangles.size();
  xvals.resize(ncells);
  yvals.resize(ncells);
  if (withBounds)
    {
      xbounds.resize(3 * ncells);
      ybounds.resize(3 * ncells);
    }

  size_t i = 0;
  for (const Triangle &t : triangles)
    {
      const auto center = circumCenterMean(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
      xvals[i] = center.longitude();
      yvals[i] = center.latitude();
      if (withBounds)
        for (size_t k = 0; k < 3; ++k)
          {
            xbounds[i * 3 + k] = vertices[t[k]].longitude();
            ybounds[i * 3 + k] = vertices[t[k]].latitude();
          }

      i++;
    }

  return ncells;
}

#ifdef TEST_ICO
#include <chrono>

static inline double
rad2deg(double v)
{
  return v * 180.0 / M_PI;
}

void
subDivdeAndVertexForEdgeTest(int subdivisions)
{
  icosahedron::init();
  auto vertices = icosahedron::vertices;
  auto triangles = icosahedron::triangles;

  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < subdivisions; ++i) triangles = subdivide(vertices, triangles);
  auto end = std::chrono::steady_clock::now();
  auto elapsedM = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

  const size_t expectedNumTriangles = std::pow(4, subdivisions) * 20;
  const size_t expectedNumVerticies = std::pow(4, subdivisions) * 10 + 2;
  if (triangles.size() != expectedNumTriangles)
    {
      std::cout << "ERROR: unexpected number of triangles for T: std::map, expected " << expectedNumTriangles << " got "
                << triangles.size() << std::endl;
    }
  if (vertices.size() != expectedNumVerticies)
    {
      std::cout << "ERROR: unexpected number of verticies for T: std::map, expected " << expectedNumVerticies << " got "
                << vertices.size() << std::endl;
    }
}

int
main(void)
{
  int numSubdivisions = 0;
  subDivdeAndVertexForEdgeTest(numSubdivisions);

  auto mesh = makeIcosphere(numSubdivisions);
  auto &vertices = mesh.first;
  auto &triangles = mesh.second;

  if (1)
    {
      for (const Vertex &v : vertices)
        {
          auto lon = v.longitude();
          auto lat = v.latitude();
          fprintf(stderr, "xyz:%g %g %g   lon:%g lat:%g\n", v[0], v[1], v[2], rad2deg(lon), rad2deg(lat));
        }
      for (Triangle &t : triangles) { fprintf(stderr, "index: %zu %zu %zu\n", t[0], t[1], t[2]); }
      for (Triangle &t : triangles)
        {
          auto center = circumCenterMean(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
          auto lon = center.longitude();
          auto lat = center.latitude();
          fprintf(stderr, "center: %g %g\n", rad2deg(lon), rad2deg(lat));
        }
      for (Triangle &t : triangles)
        {
          double lon, lat;
          printf(">\n");
          for (int i = 0; i < 3; ++i)
            {
              lon = vertices[t[i]].longitude();
              lat = vertices[t[i]].latitude();
              printf("   %g  %g\n", rad2deg(lon), rad2deg(lat));
            }
          lon = vertices[t[0]].longitude();
          lat = vertices[t[0]].latitude();
          printf("   %g  %g\n", rad2deg(lon), rad2deg(lat));
        }
    }
  fprintf(stderr, "vertices %zu\n", vertices.size());
  fprintf(stderr, "triangles %zu\n", triangles.size());
}
#endif
